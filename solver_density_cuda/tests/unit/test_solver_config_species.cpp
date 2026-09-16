// =============================================================================
// test_solver_config_species.cpp — solverConfig::read() まで通す化学種/凝縮/トレーサ設定の検査 (GPU 不要)
//   solverConfig::read() は runtime_error を捕まえて exit(EXIT_FAILURE) するため、各ケースは本プログラム自身を
//   `--read <dir>` で子プロセスとして起動し、終了コードと出力 (解決した condGasSpecies / 名前) で判定する。
//   ケース (codex 2026-09-16 result M7 / M6):
//     ok-name       condensationSpecies: H2O (species [MIXDRY, H2O]) → condGasSpecies=1, name H2O
//     ok-index      condGasSpecies: 1 だけ → name H2O (index 指定のみ)
//     ok-case       condensationSpecies: h2o (大文字小文字無視)
//     ok-both       condensationSpecies: H2O + condGasSpecies: 1 (一致)
//     mismatch      condensationSpecies: H2O + condGasSpecies: 0 → エラー
//     range         condGasSpecies: 5 → エラー (nSpecies=2)
//     unknown       condensationSpecies: XYZ → エラー
//     model-h2o0    species H2O × condModel 0 (N2 物性) → エラー
//     model-n2-1    species [N2,H2O] condensationSpecies N2 × condModel 1 → エラー
//     single        species [H2O] (nSpecies 1) × condGasSpecies 0 → エラー (carrier は 2 種以上)
//     single-pure   species [N2] pure (index -1) × condModel 0 → OK (pure N2 TP)
//     single-pure-x species [H2O] pure × condModel 0 → エラー (物質不一致)
//     tracer-dual-s0      tracer exhaust × dualTime 1 × passiveScalarScheme 0 (旧経路) → エラー (物理時間項なし)
//     tracer-dual-default tracer exhaust × dualTime 1 (既定 passiveScalarScheme 1: 受動種 BDF あり, Phase B) → OK
//     tracer-bogus  tracer: bogus → エラー
//     tracer-ok     tracer exhaust × dualTime 0 → OK
//
// ビルド/実行:
//   g++ -O1 -std=c++17 -I solver_density_cuda solver_density_cuda/tests/unit/test_solver_config_species.cpp \
//       solver_density_cuda/input/solverConfig.cpp -lyaml-cpp -o /tmp/test_solver_config_species && /tmp/test_solver_config_species
// 規約: [PASS]/[FAIL] を出し、失敗があれば非ゼロ終了。
// =============================================================================
#include "input/solverConfig.hpp"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <unistd.h>

static int g_fail = 0;
static void check(bool ok, const std::string& what)
{
    std::printf("%s %s\n", ok ? "[PASS]" : "[FAIL]", what.c_str());
    if (!ok) ++g_fail;
}

// 最小構成の solverConfig.yaml (case/44 run_0190 の config を元に species / condensation / tracer / dualTime を差し替える)
static std::string makeConfig(const std::string& species, const std::string& cond, const std::string& tracer, int dualTime,
                              const std::string& extraDeltaT = "")
{
    std::string s;
    s += "mesh: {meshFormat: \"hdf5\", discretization: \"node\", isAxisymmetric: 1, axisCentroidShift: 1, meshFileName: \"nozzle.h5\", valueFileName: \"nozzle.h5\"}\n";
    s += "gpu: 1\nsolver: \"SLAU\"\n";
    s += "physProp: {isCompressible: 1, thermalMethod: 2, viscMethod: 0, ro: 1.2, visc: 0.0, thermCond: 0.0, cp: 1220.7, gamma: 1.31526,\n";
    s += "           species: [" + species + "], speciesDBFile: \"species_db.yaml\", thermoHrefTemp: 298.15" + (tracer.empty() ? "" : ", tracer: " + tracer) + "}\n";
    s += "time:\n  unsteady: 0\n  dualTime: " + std::to_string(dualTime) + "\n  last: {control: 0, nStepOuter: 10}\n";
    s += "  deltaT: {control: 1, dt: 1e-8, cfl: 2.0, cfl_pseudo: 2.0, dt_min: 1e-9, dt_max: 0.001, blockDPLUR: 1, lowMachPrecond: 0, detectNaN: 1" + extraDeltaT + "}\n";
    s += "  outStepStart: 0\n  outStepInterval: 10\n  timeIntegration: 11\n  nStepInner: 5\n";
    s += "space: {convMethod: 1, limiter: 2}\nturbulence: {model: \"none\"}\ninitial: \"uniform_p101325_u10\"\n";
    if (!cond.empty()) s += "condensation: {" + cond + "}\n";
    s += "output: {level: 1}\n";
    return s;
}

struct Case { const char* name; std::string species; std::string cond; std::string tracer; int dualTime; bool expectOk; std::string expectOut; std::string extraDeltaT = ""; };

int main(int argc, char** argv)
{
    // ---- 子プロセスモード: 指定ディレクトリで read() し、解決結果を 1 行出す ----
    if (argc == 3 && std::string(argv[1]) == "--read") {
        if (chdir(argv[2]) != 0) { std::fprintf(stderr, "chdir failed\n"); return 2; }
        solverConfig cfg;
        cfg.read("solverConfig.yaml");   // 失敗時は内部で exit(EXIT_FAILURE)
        std::printf("RESOLVED condGasSpecies=%d name=%s tracer=%s\n", cfg.condGasSpecies, cfg.condGasSpeciesName.c_str(), cfg.tracer.c_str());
        return 0;
    }

    const std::string base = "/tmp/test_solver_config_species";
    const std::string self = argv[0];
    const std::string condH2O = "condensation: 1, nCondSpecies: 1, condModel: 1, condKantrowitz: 1, condLimiterMode: 1";
    const std::string condN2  = "condensation: 1, nCondSpecies: 1, condModel: 0, condKantrowitz: 1, condLimiterMode: 1";
    const Case cases[] = {
        {"ok-name",       "MIXDRY, H2O", condH2O + ", condensationSpecies: H2O",                     "", 0, true,  "condGasSpecies=1 name=H2O"},
        {"ok-index",      "MIXDRY, H2O", condH2O + ", condGasSpecies: 1",                            "", 0, true,  "condGasSpecies=1 name=H2O"},
        {"ok-case",       "MIXDRY, H2O", condH2O + ", condensationSpecies: h2o",                     "", 0, true,  "condGasSpecies=1 name=H2O"},
        {"ok-both",       "MIXDRY, H2O", condH2O + ", condensationSpecies: H2O, condGasSpecies: 1",  "", 0, true,  "condGasSpecies=1 name=H2O"},
        {"mismatch",      "MIXDRY, H2O", condH2O + ", condensationSpecies: H2O, condGasSpecies: 0",  "", 0, false, "disagrees"},
        {"range",         "MIXDRY, H2O", condH2O + ", condGasSpecies: 5",                            "", 0, false, "out of range"},
        {"unknown",       "MIXDRY, H2O", condH2O + ", condensationSpecies: XYZ",                     "", 0, false, "is not in physProp.species"},
        {"model-h2o0",    "MIXDRY, H2O", condN2  + ", condensationSpecies: H2O",                     "", 0, false, "does not match condModel 0"},
        {"model-n2-1",    "N2, H2O",     condH2O + ", condensationSpecies: N2",                      "", 0, false, "does not match condModel 1"},
        {"single",        "H2O",         condH2O + ", condGasSpecies: 0",                            "", 0, false, "requires >=2 species"},
        {"single-pure",   "N2",          condN2,                                                     "", 0, true,  "condGasSpecies=-1 name=N2"},
        {"single-pure-x", "H2O",         condN2,                                                     "", 0, false, "does not match condModel 0"},
        {"tracer-dual-s0",      "MIXDRY, H2O", "",                                                   "exhaust", 1, false, "requires passiveScalarScheme 1", ", passiveScalarScheme: 0"},
        {"tracer-dual-default", "MIXDRY, H2O", "",                                                   "exhaust", 1, true,  "tracer=exhaust"},
        {"tracer-bogus",  "MIXDRY, H2O", "",                                                         "bogus",   0, false, "must be 'none' or 'exhaust'"},
        {"tracer-ok",     "MIXDRY, H2O", "",                                                         "exhaust", 0, true,  "tracer=exhaust"},
    };
    for (const Case& c : cases) {
        const std::string dir = base + "/" + c.name;
        std::system(("mkdir -p " + dir).c_str());
        { std::ofstream f(dir + "/solverConfig.yaml"); f << makeConfig(c.species, c.cond, c.tracer, c.dualTime, c.extraDeltaT); }
        const std::string out = dir + "/out.txt";
        const int rc = std::system((self + " --read " + dir + " > " + out + " 2>&1").c_str());
        const bool ok = (rc == 0);
        std::string text, line;
        { std::ifstream f(out); while (std::getline(f, line)) text += line + "\n"; }
        const bool found = text.find(c.expectOut) != std::string::npos;
        check(ok == c.expectOk && found,
              std::string(c.name) + ": expected " + (c.expectOk ? "OK" : "ERROR") + " containing '" + c.expectOut + "' -> "
              + (ok ? "OK" : "ERROR") + (found ? " (found)" : " (NOT found)"));
        if (!found) {
            // 判定に使った行を表示 (Configuration Error / RESOLVED)
            std::string::size_type p = text.find("Configuration Error");
            if (p == std::string::npos) p = text.find("RESOLVED");
            if (p != std::string::npos) std::printf("       %s\n", text.substr(p, text.find('\n', p) - p).c_str());
        }
    }
    std::printf("%s (%d failures)\n", g_fail == 0 ? "ALL PASS" : "FAILED", g_fail);
    return g_fail == 0 ? 0 : 1;
}
