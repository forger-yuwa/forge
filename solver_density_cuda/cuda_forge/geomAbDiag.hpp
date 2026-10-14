#pragma once

// V0 の評価の経路 (plans/active/architecture-float-state-double-geometry.md §4.2c、段 ①)。既定 off。
//
// 状態を固定して、粘性 (viscousFlux_d) と SST k/ω の拡散 (scalar_diffusion_multi_d) の面の流束を、
// 座標の差の取り方だけを変えて比べるための経路。main の runGeomAbDiag が組立 (assembleResidualPre + Post) を
// 1 回だけ通し (= 1 step 目の組立)、本番の起動の直前・直後にここのフックが呼ばれる。境界条件・EOS・勾配・
// 乱流モデル・commit は腕ごとにやり直さない。h5 を書いたら時間更新へ進まず終了する。
//
//   FORGE_DIAG_GEOMAB_DUMP=<h5> (float のビルド)
//     - 本番の起動の直前の入力 (カーネルが読むセル配列・面の幾何・接続・float の座標・e32) を記録する。
//     - 同じ入力で、旧腕 (今の座標の差 ccx[ic1] − ccx[ic0]) と新腕 (e32 = double の座標の差を 1 回丸めた値) を、
//       別の残差バッファに向けて 1 回ずつ評価し、atomicAdd の直前の面の流束を面の番号つきで書く。
//     - 旧腕の残差は「本番の起動の直前の残差の写し」から積み、本番の起動の直後の残差とビット単位で照合する
//       (引き算の丸めを挟まずに「寄与」を比べるため)。atomicAdd の順序で一致しない分の目安として、旧腕をもう
//       1 回同じ条件で評価した差も数える。
//   FORGE_DIAG_GEOMAB_REF=<h5> (FP64 のビルド)
//     - <h5> は DUMP の出力。同じ位置で、記録した入力 (float を double に広げたもの) でセル配列と面の幾何
//       (fx・sx..ss) を上書きし、FP64 のビルド自身の e (double の座標の差) で同じカーネルを 1 回評価する。
//     - 面の流束を <h5 の末尾 .h5 を _ref.h5 に替えた名前> に書く (DUMP の出力は書き換えない)。
//
// 面の流束の並び: スカラー拡散は [n, nPlanes] (s = 0..n-1 は k, ω の順)、粘性は [6, nPlanes]
// (運動量 x, y, z・エネルギー・その内の熱・仕事。viscousFlux_d の引数の説明)。評価しない面は NaN。

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"
#include "cuda_forge/cudaConfig.cuh"

struct ScalarTransportDesc;

namespace geomAbDiag {

// 腕の起動 (本番と同じ引数の組み立てで、書き先・e・面の流束の出力先だけを差し替える)。
// フックの呼び出し側 (viscousFlux_d.cu / scalarTransport_d.cu) が自分の起動関数を渡す。ここから直接参照しないのは、
// 変換器 (convertGmshToForge) が libcuda_forge から scalarTransport_d.o を取り込むとき、粘性のオブジェクト経由で
// conjugateWall (forge にしか無い) まで芋づるに引き込まないため。
using ViscousArmFn = void (*)(solverConfig&, cudaConfig&, mesh&, variables&, flow_float* const*,
                              const flow_float*, const flow_float*, const flow_float*, flow_float*);
using ScalarArmFn  = void (*)(solverConfig&, cudaConfig&, mesh&, variables&, const ScalarTransportDesc*, int,
                              flow_float* const*, flow_float* const*,
                              const flow_float*, const flow_float*, const flow_float*, flow_float*);

// 0 = off、1 = DUMP、2 = REF。環境変数を初回に読む (両方あれば起動を止める)。
int mode();

// 組立を 1 回通す間だけ立つ (main の runGeomAbDiag が立てて下ろす)。下りている間、フックは呼ばれない。
bool armed();
void arm(bool on);

// 開始時に 1 回: 接続の記録、e32 と今の差の照合、REF では DUMP の出力との大きさ・接続の一致を確かめる。
// 不一致なら理由を出して false。
bool begin(solverConfig& cfg, mesh& msh, variables& var);

// SST k/ω の融合拡散 (scalarTransportResidualMulti_d) の本番の起動の直前・直後。armFn = scalarDiffusionMultiArm_d
void scalarBefore(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var,
                  const ScalarTransportDesc* descs, int n, ScalarArmFn armFn);
void scalarAfter(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var,
                 const ScalarTransportDesc* descs, int n);

// 内部面の粘性流束 (viscousFlux_d_wrapper) の本番の起動の直前・直後。armFn = viscousFluxInternalArm_d
void viscousBefore(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, ViscousArmFn armFn);
void viscousAfter(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var);

// 記録を h5 に書く。戻り値は main の終了コード。
int finish(solverConfig& cfg, mesh& msh, variables& var);

}  // namespace geomAbDiag
