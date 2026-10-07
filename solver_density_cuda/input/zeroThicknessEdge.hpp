#pragma once
// =============================================================================
// 厚さ 0 の板の自由端の近傍で速度の再構成を節点値へ寄せる (space.zeroThicknessEdgeVelocity) の起動時処理。
// plans/active/convection-zero-thickness-edge-reconstruction.md §4。
//
// 端の判定 (端集合 E・距離 rings 以内の集合 S_2) は**ソルバでなく前処理**
// (tools/mark_zero_thickness_edges.py) が行い、メッシュ h5 の /AUX/w_recon_vel に節点 (node の CV) ごとの重み
// w∈[0,1] を書く。ソルバはそれを読んで SLAU の内部面で u_f = u_i + w_i ψ∇u_i·r_if にするだけ。
//
// 置き場所を /VALUE でなく /AUX にする理由: restart_field.py / interp_field.py は /VALUE だけを書き換え、
// 他のグループは宛先 (= 新しい入力 h5) のものを残す。w はメッシュの位相とタグだけで決まる量なので、
// 同じ格子の継続ではそのまま残り、格子を変えたら新しい h5 に /AUX が無い → 起動時エラーで作り直しを強制できる。
// =============================================================================
#include <string>
#include <vector>
#include "flowFormat.hpp"
#include "input/solverConfig.hpp"
#include "mesh/mesh.hpp"

// 起動時に 1 回 (readBcondConfig の後、bcond 種別が揃ってから) 呼ぶ。
//   (1) 起動エコーを 1 行出す (無効でも出す。stage_manifest・来歴の正本)。
//   (2) 有効なら契約を検査: node・solver SLAU/SLAU2・gpu 1・周期境界なし・非軸対称。
//   (3) メッシュ h5 (cfg.meshFileName) の /AUX/w_recon_vel を読み、型 (float32)・長さ (= nCells)・値域 [0,1]・有限性、
//       属性 field_sha256 (値の SHA-256) と mesh_signature (読み込んだ格子の座標・内部面・境界面の接続から再計算する
//       格子署名、版 zte-mesh-sig-v1。定義は .cpp の meshSignature と道具の mesh_signature()) を照合する。
//   (4) 起動ログ: w<1 の節点数・w の min/max・座標範囲・対象体積・対象の内部面数・生成元の属性。
//       env FORGE_EDGE_MASK_PROBE="<id>,<id>,..." で指定節点の w を出す。
// 違反は理由を出して exit(1)。戻り値: 有効なら w [nCells] (device へ上げる)、無効なら空。
// 起動記録 (main.cpp の forge_launches.jsonl) に書く照合済みの値。ハッシュは属性の転記でなく再計算した値。
struct ZteLaunchInfo {
    int enabled = 0;
    std::string field;                  // "w_recon_vel"
    std::string meshSignature;          // 読み込んだ格子から再計算した格子署名 (16 進 64 桁)
    std::string meshSignatureVersion;   // "zte-mesh-sig-v1"
    std::string fieldSha256;            // 読み込んだ w の SHA-256 (float32 のバイト列)
    long long nZero = 0, nLt1 = 0, nNodes = 0;
};

std::vector<flow_float> zteLoadVelocityWeight(const solverConfig& cfg, const mesh& msh, ZteLaunchInfo* info = nullptr);
