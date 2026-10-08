#pragma once
// 毎更新の EOS 床事象のカウンタ (出力専用; config `output: {floorEvents: 1}`、既定 0 = 何もしない)。
// 定義の正本は plans/active/tooling-sern-te-wake-grid.md §4「床の判定」(2026-10-08 のやり直し)。
//
// 数えるのは **EOS の床を必要とした事象** (定常陰解法では床は作業配列だけに入り、commit は床前の基準
// roN/roeN + dq なので、「保存量に残った補正」は構造上 0 件になる)。床の種類別の述語で数え、補正量は付帯情報:
//   温度: EOS と同じ密度 (床後 ρ)・組成での e_in < e_mix(T_min)。EOS の最終温度が下限かどうかとの食い違いは別列。
//   密度: ρ_in < roMin。補正量 Δρ = roMin − ρ_in。
//   圧力: クランプ直前の P_raw < pMin。補正量 ΔP = pMin − P_raw (TP では roe に触れないので直接の Δ(ρE) は 0)。
// 補正量による足切りはしない (除外閾値 0)。述語が成立した事象だけを数えるので、通常状態の roe 再構成の丸めは数えない。
// 判定はすべて「境界ピン (no-slip・軸・等温壁) の後・EOS の床の前」の状態 = dependentVariables_d の入力で行う。
// 記録形式・時点の契約は floorEvents_d.cu 冒頭と procedures/solver-settings.md の「output」節。

#include <string>

#include "flowFormat.hpp"
#include "cuda_forge/cudaConfig.cuh"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

constexpr int FLOOR_EVENT_IDS = 8;   // 種類ごとに記録する節点 ID の上限 (超えたら件数は数え続け、ID 列だけ切り詰め印)

// ID を記録する種類 (FloorEventAcc::ids の添字)
enum FloorEventIdKind { FE_ID_T = 0, FE_ID_RHO = 1, FE_ID_P = 2, FE_ID_NEAR = 3, FE_ID_KINDS = 4 };

// 1 回の EOS 呼び出しの集計 (device)。呼び出しごとに 0 へ戻す。[0] = 実節点 (ic < nCells)、[1] = ghost。
struct FloorEventAcc {
    unsigned long long nNonfinite[2];   // 入力の保存量 (ro, roU, roe) が非有限 (床の述語から外す)
    unsigned long long nT[2];           // 温度床
    unsigned long long nRho[2];         // 密度床
    unsigned long long nP[2];           // 圧力床
    unsigned long long nTUneval;        // 実節点: 温度の述語を評価しない経路 (二相・EOS 拘束形平衡・二相反転の失敗)
    unsigned long long nTMismatch;      // 実節点: 述語と「EOS の最終温度が下限」の食い違い
    unsigned long long nNear;           // 実節点: EOS の最終温度 ≤ T_min + nearBand (補助。床事象の節点も含む)
    double sum_dRhoE_T;                 // 実節点: Σ ρ (e_mix(T_min) − e_in)   [J/m³]
    double sum_dRho;                    // 実節点: Σ (roMin − ρ_in)            [kg/m³]
    double sum_dP;                      // 実節点: Σ (pMin − P_raw)            [Pa]
    unsigned long long max_dRhoE_T;     // 非負 double のビット列 (atomicMax で最大を取る)
    unsigned long long max_dRho;
    unsigned long long max_dP;
    unsigned int nIds[FE_ID_KINDS];     // ID の記録を試みた数 (> FLOOR_EVENT_IDS なら切り詰め)
    long long ids[FE_ID_KINDS][FLOOR_EVENT_IDS];
};

// カーネル引数。acc == nullptr なら計数しない (無効時の EOS は従来と同じ演算列)。
struct FloorEventDev {
    FloorEventAcc* acc;
    const double* hTminMass;   // 化学種ごとの h_s(T_min)/MW_s [J/kg] (thermo_cph_mix と同じ項)。TP のみ
    double tFloorTP;           // TP の温度下限 (DEPVAR_TMIN)
    double nearBand;           // 床近傍の幅 [K]
};

#if defined(__CUDACC__)
__device__ inline void feAtomicMaxNonneg(unsigned long long* m, double v)
{
    if (v > 0.0) atomicMax(m, (unsigned long long)__double_as_longlong(v));
}

__device__ inline void feRecordId(FloorEventAcc* a, int kind, geom_int ic)
{
    const unsigned int s = atomicAdd(&a->nIds[kind], 1u);
    if (s < (unsigned int)FLOOR_EVENT_IDS) a->ids[kind][s] = (long long)ic;
}

// 入力の保存量の有限性と密度床。戻り値: 入力が有限か (偽なら以降の述語を評価しない)。
__device__ inline bool feInputAndDensity(const FloorEventDev& fe, geom_int ic, bool real,
                                         flow_float ro_in, flow_float rux, flow_float ruy, flow_float ruz, flow_float roe_in,
                                         flow_float roMin)
{
    FloorEventAcc* a = fe.acc;
    const int w = real ? 0 : 1;
    if (!(isfinite(ro_in) && isfinite(rux) && isfinite(ruy) && isfinite(ruz) && isfinite(roe_in))) {
        atomicAdd(&a->nNonfinite[w], 1ull);
        return false;
    }
    if (ro_in < roMin) {
        atomicAdd(&a->nRho[w], 1ull);
        if (real) {
            const double d = (double)roMin - (double)ro_in;
            atomicAdd(&a->sum_dRho, d); feAtomicMaxNonneg(&a->max_dRho, d);
            feRecordId(a, FE_ID_RHO, ic);
        }
    }
    return true;
}

// 圧力床: クランプ直前の P_raw < pMin。
__device__ inline void fePressure(const FloorEventDev& fe, geom_int ic, bool real, double P_raw, double pMin)
{
    if (!(P_raw < pMin)) return;
    FloorEventAcc* a = fe.acc;
    atomicAdd(&a->nP[real ? 0 : 1], 1ull);
    if (real) {
        const double d = pMin - P_raw;
        atomicAdd(&a->sum_dP, d); feAtomicMaxNonneg(&a->max_dP, d);
        feRecordId(a, FE_ID_P, ic);
    }
}

// 温度床: 述語 pred (e_in < e_mix(T_min)) と EOS の最終温度。dRhoE = ρ (e_mix(T_min) − e_in)。
__device__ inline void feTemperature(const FloorEventDev& fe, geom_int ic, bool real, bool pred,
                                     bool finalAtFloor, double T_final, double T_floor, double dRhoE)
{
    FloorEventAcc* a = fe.acc;
    if (pred) {
        atomicAdd(&a->nT[real ? 0 : 1], 1ull);
        if (real) {
            atomicAdd(&a->sum_dRhoE_T, dRhoE); feAtomicMaxNonneg(&a->max_dRhoE_T, dRhoE);
            feRecordId(a, FE_ID_T, ic);
        }
    }
    if (!real) return;
    if (pred != finalAtFloor) atomicAdd(&a->nTMismatch, 1ull);
    if (T_final <= T_floor + fe.nearBand) {
        atomicAdd(&a->nNear, 1ull);
        feRecordId(a, FE_ID_NEAR, ic);
    }
}

__device__ inline void feTemperatureUnevaluated(const FloorEventDev& fe, bool real)
{
    if (real) atomicAdd(&fe.acc->nTUneval, 1ull);
}
#endif

// ---- host 側 (floorEvents_d.cu) ------------------------------------------------------------------
// 起動時: 有効なら記録ファイル floor_events.csv を開き (追記)、セッション開始行を書く。cfgFnv は solverConfig.yaml の FNV-1a。
void floorEventsOpen(const solverConfig& cfg, const mesh& msh, const std::string& cfgFnv);
bool floorEventsEnabled();
// 計数の文脈: 初期化の EOS / 時間ステップ内の EOS / それ以外 (診断の再組立など) / 終了時の監査
enum class FloorEventPhase { Init, Step, Aux, Audit };
void floorEventsSetPhase(FloorEventPhase phase, int step);
// EOS カーネルの引数 (無効なら acc == nullptr)。呼ぶたびに集計を 0 へ戻す。
FloorEventDev floorEventsKernelArgs();
// EOS カーネルの後: 集計を読み、1 行書く (無効なら何もしない)。
void floorEventsAfterEos();
// 終了時の監査: 最後の更新結果 Q_N を、コピー上で同じ前処理 (境界ピン) と同じ EOS (同じ判定) にかける。
// 通常の配列には触れない (var.c_d の該当ポインタを一時的にコピーへ差し替える)。
void floorEventsAudit(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int lastStep);
// セッション終了行を書いて閉じる。
void floorEventsClose(int lastStep);
// 試験用 (既定 off、計測の閉じ方の小型 A/B: plan §5.1 #2)。step の更新の後で 1 節点の保存量だけを書き換える。
//   FORGE_FLOOR_TEST_INJECT=step:node:de:<J/kg>       roe だけを ρ (e_mix(T_min) + ek + de) に (密度・運動量・組成は不変)
//   FORGE_FLOOR_TEST_INJECT=step:node:rhoscale:<倍率>  その節点の保存量を全部同じ倍率に (速度・e・組成は不変、ρ と p が下がる)
// FORGE_FLOOR_AUDIT_VERIFY=1 で監査の前後の全 cell 配列をバイト比較して表示する。
void floorEventsTestInject(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int step);
