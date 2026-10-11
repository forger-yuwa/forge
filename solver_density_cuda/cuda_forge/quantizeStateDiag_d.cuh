#pragma once

// 診断 FORGE_DIAG_QUANTIZE_STATE の 1 値の処理と、自己点検の集計・書式
// (plans/active/architecture-float-state-double-geometry.md §5.1 #19・§6.36、既定 off)。
// 制御 (環境変数・拒否・フックの時点) は quantizeStateDiag.hpp。
//
// 単体試験 (tests/unit/test_quantize_state.cu) からも読むので、依存は flowFormat.hpp と std だけにする。
//   - quantizeStateValue        : カーネル (quantizeStateDiag_d.cu) が 1 CV の 1 量に当てる処理そのもの。
//   - quantizeStateRepresentable: 自己点検の判定 ((double)(float)x == x)。
//   - quantizeStateFormat       : 自己点検の 1 行の本文 (量ごとの個数の並び)。
//   - quantizeStateParseOnly    : FORGE_DIAG_QUANTIZE_STATE_ONLY の値の検査 (§5.1 #20)。
//   - quantizeStateCountChanged : 丸めの直前と直後の生のビットの比較 (選ばなかった量が変わらないことの確認、§5.1 #20)。

#include "flowFormat.hpp"

#include <cstddef>
#include <cstring>
#include <string>
#include <vector>

#ifdef __CUDACC__
  #define QUANTIZE_STATE_HD __host__ __device__ inline
#else
  #define QUANTIZE_STATE_HD inline
#endif

// 1 値を float32 に丸めて flow_float に戻す。丸めは最近接偶数 (C の変換。device は cvt.rn.f32.f64 で、-ftz を付けない
// ビルドなので float32 の非正規化数も保つ)。FP64 のビルドでは (double)(float)x。float のビルドでは恒等 (このビルドは
// 起動時に拒否するので呼ばれない)。float32 で表せる値はそのまま返す (2 回当てても変わらない)。
QUANTIZE_STATE_HD flow_float quantizeStateValue(flow_float x)
{
    return static_cast<flow_float>(static_cast<float>(x));
}

// float32 で表せるか ((double)(float)x == x)。NaN は false (表せないに数える)。±0 は true。
inline bool quantizeStateRepresentable(double x)
{
    return static_cast<double>(static_cast<float>(x)) == x;
}

// [i0, i1) の表せない値の個数。skip があれば (*skip)[i] == 1 の点を数えない (node の入口のピン scalarDirichletPin 用)。
template <typename T, typename P = T>
inline long long quantizeStateCountNotRepresentable(const std::vector<T>& v, size_t i0, size_t i1,
                                                     const std::vector<P>* skip = nullptr)
{
    long long n = 0;
    if (i1 > v.size()) i1 = v.size();
    if (skip != nullptr && i1 > skip->size()) i1 = skip->size();
    for (size_t i = i0; i < i1; ++i) {
        if (skip != nullptr && (*skip)[i] == static_cast<P>(1)) continue;
        if (!quantizeStateRepresentable(static_cast<double>(v[i]))) ++n;
    }
    return n;
}

// 自己点検の 1 群 (Q・Q_N・Q_M・roY_N・roY_M)。names と counts は同じ長さ。
struct QuantizeStateCountGroup {
    std::string label;                 // 例 "Q"、"Q_N"
    std::vector<std::string> names;    // 例 "ro"、"roUx"
    std::vector<long long> counts;     // 量ごとの表せない値の個数
};

// "ro 0 roUx 0 ... (Q) / roN 0 ... (Q_N) / 合計 0" を作る。total に合計を返す (nullptr なら返さない)。
inline std::string quantizeStateFormat(const std::vector<QuantizeStateCountGroup>& groups, long long* total)
{
    std::string s;
    long long sum = 0;
    for (const auto& g : groups) {
        for (size_t k = 0; k < g.names.size(); ++k) {
            s += g.names[k] + " " + std::to_string(g.counts[k]) + " ";
            sum += g.counts[k];
        }
        s += "(" + g.label + ") / ";
    }
    s += "合計 " + std::to_string(sum);
    if (total != nullptr) *total = sum;
    return s;
}

// 時点の説明 (行の末尾の括弧に入れる)。
inline std::string quantizeStateWhenOuter() { return "updateVariablesOuter の直後"; }
inline std::string quantizeStateWhenSpecies(int stepNow)
{
    return "step " + std::to_string(stepNow) + " の化学種の更新の頭 (speciesUpdateOuter の直後)";
}

// 自己点検の行 (改行なし)。step は点検の対象の step (iStep + 1)。
//   "[quantizeState] step N: 表せない値 <body> (全 CV <nAll>、<when>)"
inline std::string quantizeStateCheckLine(int step, const std::string& body, long long nAll, const std::string& when)
{
    return "[quantizeState] step " + std::to_string(step) + ": 表せない値 " + body + " (全 CV " + std::to_string(nAll) +
           "、" + when + ")";
}
// 合計が 0 でないときに続けて出す内点の行 (改行なし)。nPin >= 0 なら、入口のピンの節点を除いた内点の個数を続ける
// (組立の中で状態を書く入口の Dirichlet のピンと、丸めの漏れを見分けるため)。
//   "[quantizeState] step N: 注意: 表せない値がある (0 を期待)。内点 (ic < nCells <nInt>) に限ると <bodyInt>
//    。入口のピン (scalarDirichletPin 1、<nPin> 節点) を除く内点では <bodyNoPin>"
inline std::string quantizeStateInteriorLine(int step, long long nInt, const std::string& bodyInt,
                                             long long nPin = -1, const std::string& bodyNoPin = "")
{
    std::string s = "[quantizeState] step " + std::to_string(step) + ": 注意: 表せない値がある (0 を期待)。内点 (ic < nCells " +
                    std::to_string(nInt) + ") に限ると " + bodyInt;
    if (nPin >= 0)
        s += "。入口のピン (scalarDirichletPin 1、" + std::to_string(nPin) + " 節点) を除く内点では " + bodyNoPin;
    return s;
}

// ---- FORGE_DIAG_QUANTIZE_STATE_ONLY (丸める量を選ぶ、§5.1 #20・§6.40) ----

// 値 ("roe"、"roe,roK") を読む。all は丸めの対象の全量 (起動時の構成で決まる。ro..roe、SST なら roK・roOmega、化学種なら
// roY0..)。成功なら true を返し、*selected を all と同じ長さにして選んだ量を true にする (並びは all のまま。値の中の順は
// 問わない)。失敗なら false を返し、*err に理由を入れる。拒否するもの (どれも黙って違う試験にしないため):
//   - 空の値 (台本の変数の空展開で、黙って全量の丸めにならないように)
//   - 空の要素 ("roe,"、",roe"、"roe,,ro")
//   - 対象に無い名前 (この構成に無い量。roK を SST でない run で選ぶ等)。前後の空白は取り除かない (" roe" は対象に無い名前)
//   - 同じ名前の重複 ("roe,roe")
inline bool quantizeStateParseOnly(const std::string& s, const std::vector<std::string>& all,
                                   std::vector<bool>* selected, std::string* err)
{
    std::vector<bool> sel(all.size(), false);
    if (s.empty()) {
        if (err != nullptr) *err = "値が空";
        return false;
    }
    size_t pos = 0;
    int k = 0;
    while (true) {
        const size_t comma = s.find(',', pos);
        const std::string tok = s.substr(pos, (comma == std::string::npos) ? std::string::npos : comma - pos);
        ++k;
        if (tok.empty()) {
            if (err != nullptr) *err = std::to_string(k) + " 番目の要素が空";
            return false;
        }
        size_t idx = all.size();
        for (size_t i = 0; i < all.size(); ++i) if (all[i] == tok) { idx = i; break; }
        if (idx == all.size()) {
            if (err != nullptr) *err = "'" + tok + "' は丸めの対象に無い";
            return false;
        }
        if (sel[idx]) {
            if (err != nullptr) *err = "'" + tok + "' が重複している";
            return false;
        }
        sel[idx] = true;
        if (comma == std::string::npos) break;
        pos = comma + 1;
    }
    if (selected != nullptr) *selected = sel;
    return true;
}

// 起動時の行の「丸める量」の欄 (選んだとき)。
//   "選んだ量だけ (roe)、FORGE_DIAG_QUANTIZE_STATE_ONLY。丸めない対象: ro, roUx, roUy, roUz"
inline std::string quantizeStateOnlyText(const std::vector<std::string>& all, const std::vector<bool>& selected)
{
    std::string on, off;
    for (size_t i = 0; i < all.size(); ++i) {
        std::string& t = (i < selected.size() && selected[i]) ? on : off;
        t += (t.empty() ? "" : ", ") + all[i];
    }
    return "選んだ量だけ (" + on + ")、FORGE_DIAG_QUANTIZE_STATE_ONLY。丸めない対象: " + (off.empty() ? std::string("なし") : off);
}

// a と b で生のビットが違う値の個数 (丸めのカーネルの直前と直後の比較)。NaN も同じビットなら変わっていない、+0 と −0 は
// 変わった、に数える。長さが違えば、長い方の余りを変わったに数える。
template <typename T>
inline long long quantizeStateCountChanged(const std::vector<T>& a, const std::vector<T>& b)
{
    const size_t n = (a.size() < b.size()) ? a.size() : b.size();
    long long c = (long long)((a.size() > b.size()) ? a.size() - b.size() : b.size() - a.size());
    for (size_t i = 0; i < n; ++i) c += (std::memcmp(&a[i], &b[i], sizeof(T)) != 0);
    return c;
}

// 自己点検の step で、丸めのカーネルの直前と直後に選ばなかった対象を比べた行 (改行なし)。names と counts は同じ長さ。
// total に合計を返す (nullptr なら返さない)。
//   "[quantizeState] step N: 選ばなかった量の丸め前後の変化 ro 0 roUx 0 ... / 合計 0
//    (0 が期待。全 CV <nAll>、丸めのカーネルの直前と直後の生のビットの比較)"
//   選ばなかった量が無いときは "... 選ばなかった量の丸め前後の変化 なし (選んだ量が対象のすべて)"
inline std::string quantizeStateUnselectedLine(int step, const std::vector<std::string>& names,
                                               const std::vector<long long>& counts, long long nAll, long long* total)
{
    std::string s = "[quantizeState] step " + std::to_string(step) + ": 選ばなかった量の丸め前後の変化 ";
    long long sum = 0;
    if (names.empty()) {
        s += "なし (選んだ量が対象のすべて)";
    } else {
        for (size_t k = 0; k < names.size(); ++k) {
            const long long c = (k < counts.size()) ? counts[k] : -1;
            s += names[k] + " " + std::to_string(c) + " ";
            sum += (c > 0) ? c : 0;
        }
        s += "/ 合計 " + std::to_string(sum) + " (0 が期待。全 CV " + std::to_string(nAll) +
             "、丸めのカーネルの直前と直後の生のビットの比較)";
    }
    if (total != nullptr) *total = sum;
    return s;
}

// カーネルのラッパ (quantizeStateDiag_d.cu)。q[0..n) に quantizeStateValue を当てる。既定のストリーム。
struct cudaConfig;
void quantizeStateToFloat32_d_wrapper(cudaConfig& cuda_cfg, geom_int n, flow_float* q);
