#pragma once
// condSonicModel (凝縮セルの二相 frozen 音速) の自動解決 (plans/active/condensation-kantrowitz-gamma-twophase-sonic.md §4.2)。
// 自動 (requested<0) では **検証済み構成** (TP carrier H2O, condEquilibrium 0, 全 bcond が inlet_Pressure/outflow/wall/slip/periodic)
// のときだけ 1 を返す。それ以外 (pure N2 / CPG / 平衡 / outlet_statPress・wall_isothermal 等の ghost 再構築が全蒸気 EOS の境界) は 0。
// 明示指定 (0/1) は尊重するが、未検証構成で 1 なら warn に理由を書く。依存無しの純関数 (tests/unit/test_cond_sonic.cpp で回帰試験)。
#include <string>
#include <vector>

inline int resolveCondSonicModel(int requested, int condensation, int thermalMethod, int condGasSpecies,
                                 int condModel, int condEquilibrium, const std::vector<std::string>& bcKinds,
                                 std::string& reason, std::string& warn)
{
    static const char* verified[] = {"inlet_Pressure", "outflow", "wall", "slip", "periodic"};
    std::string bad;
    for (const auto& k : bcKinds) {
        bool ok = false;
        for (const char* v : verified) if (k == v) ok = true;
        if (!ok) bad += (bad.empty() ? "" : ",") + k;
    }
    std::string why;
    if (condensation != 1)      why = "condensation off";
    else if (thermalMethod != 2) why = "thermalMethod!=2 (CPG 分岐は旧式のみ)";
    else if (condGasSpecies < 0) why = "pure-condensible (未検証)";
    else if (condModel != 1)     why = "condModel!=1 (H2O 以外は未検証)";
    else if (condEquilibrium != 0) why = "condEquilibrium!=0 (緩和形/EOS 拘束形は未検証)";
    else if (!bad.empty())       why = "未検証境界 (" + bad + ") の ghost/ピンが全蒸気 EOS";
    const bool verifiedCfg = why.empty();
    warn.clear();
    if (requested >= 0) {
        reason = "explicit condSonicModel=" + std::to_string(requested);
        if (requested == 1 && !verifiedCfg) warn = "condSonicModel=1 は未検証構成: " + why;
        return requested;
    }
    if (condensation != 1) { reason = "auto: " + why; return 0; }
    reason = verifiedCfg ? "auto: TP carrier H2O, condEquilibrium 0, 検証済み境界のみ -> 1" : ("auto: " + why + " -> 0");
    return verifiedCfg ? 1 : 0;
}
