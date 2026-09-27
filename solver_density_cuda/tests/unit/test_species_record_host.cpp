// =============================================================================
// test_species_record_host.cpp — 化学種の解決済み記録・互換性ハッシュ・入力場照合 (input/speciesDB.cpp) の単体試験 (GPU 不要)
//   plans/active/thermophysics-solver-owned-species-db.md §5.1 #3a / §6 V1
//   (1) SHA-256 の既知ベクトル (FIPS 180-4: "", "abc", 2 ブロック)
//   (2) 互換性ハッシュ: source だけ違う (同一係数を外部 DB で与える) → 一致 / N2 nasa9_low[2] +0.001 → 不一致 /
//       thermoHrefTemp 違い → 不一致 / 種の順序違い → 不一致
//   (3) 記録の書き込み: 同一内容は再利用 (上書きしない)、来歴違いは別名、記録の差分に該当係数が出る
//   (4) 入力場の照合: 一致 → 通す / 属性なし → 照合不能で拒否・env 許可で通して未検証印 / 不一致 → 該当係数を示して拒否 (env でも通さない) /
//       記録の取り違え (別 run の記録を置く) → 完全性不一致を明示 / 記録なし → 「特定不能」を明示 / 入力の未検証印は継承
//
// ビルド/実行:
//   g++ -O1 -std=c++17 -I solver_density_cuda solver_density_cuda/tests/unit/test_species_record_host.cpp
//       solver_density_cuda/input/speciesDB.cpp -lyaml-cpp -o /tmp/test_species_record_host && /tmp/test_species_record_host
//   (上の 2 行を 1 行につなげて実行)
// 規約: [PASS]/[FAIL] を出し、失敗があれば非ゼロ終了。
// =============================================================================
#include "input/speciesDB.hpp"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

static int g_fail = 0;
static void check(bool ok, const std::string& what)
{
    std::printf("%s %s\n", ok ? "[PASS]" : "[FAIL]", what.c_str());
    if (!ok) ++g_fail;
}
static bool has(const std::string& s, const std::string& sub) { return s.find(sub) != std::string::npos; }

static std::string readAll(const fs::path& p)
{
    std::ifstream f(p, std::ios::binary);
    std::ostringstream ss; ss << f.rdbuf();
    return ss.str();
}

// 内蔵 N2 を外部 DB 形式で書く (low[2] に dlow2 を足す)
static void writeN2Db(const fs::path& p, double dlow2)
{
    const auto b = speciesDB_builtin().at("N2");
    std::ofstream f(p);
    char buf[64];
    f << "N2:\n";
    std::snprintf(buf, sizeof(buf), "%.17g", b.MW); f << "  MW: " << buf << "\n";
    std::snprintf(buf, sizeof(buf), "%.17g", b.sigma_LJ); f << "  LJ_sigma: " << buf << "\n";
    std::snprintf(buf, sizeof(buf), "%.17g", b.eps_kB); f << "  LJ_eps_kB: " << buf << "\n";
    f << "  Tlo: 200.0\n  Tmid: 1000.0\n  Thi: 6000.0\n";
    f << "  nasa9_low: [";
    for (int k = 0; k < 9; ++k) { std::snprintf(buf, sizeof(buf), "%.17g", b.low[k] + (k == 2 ? dlow2 : 0.0)); f << (k ? ", " : "") << buf; }
    f << "]\n  nasa9_high: [";
    for (int k = 0; k < 9; ++k) { std::snprintf(buf, sizeof(buf), "%.17g", b.high[k]); f << (k ? ", " : "") << buf; }
    f << "]\n";
}

int main()
{
    // ---- (1) SHA-256 ----
    check(speciesDB_sha256Hex("") == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855", "sha256('')");
    check(speciesDB_sha256Hex("abc") == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad", "sha256('abc')");
    check(speciesDB_sha256Hex("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq")
          == "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1", "sha256(2 blocks)");
    {
        const std::string m(1000000, 'a');
        check(speciesDB_sha256Hex(m) == "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0", "sha256(1e6 x 'a')");
    }

    const fs::path tmp = fs::temp_directory_path() / "forge_test_species_record";
    fs::remove_all(tmp);
    fs::create_directories(tmp);
    const fs::path dbSame = tmp / "db_same.yaml", dbMod = tmp / "db_mod.yaml";
    writeN2Db(dbSame, 0.0);
    writeN2Db(dbMod, 0.001);

    // ---- (2) 互換性ハッシュ ----
    const ResolvedSpeciesDB dbB = speciesDB_resolve({"N2", "H2O"}, "");
    const ResolvedSpeciesDB dbF = speciesDB_resolve({"N2", "H2O"}, dbSame.string());
    const ResolvedSpeciesDB dbM = speciesDB_resolve({"N2", "H2O"}, dbMod.string());
    const ResolvedSpeciesDB dbR = speciesDB_resolve({"H2O", "N2"}, "");
    check(dbF.source[0] == "file" && dbB.source[0] == "builtin", "source: builtin vs file");
    const std::string hB = speciesDB_compatHash(dbB, 298.15);
    check(hB.size() == 64, "compat hash length 64");
    check(speciesDB_compatHash(dbF, 298.15) == hB, "same coefficients from external DB (source differs) -> same compat hash");
    check(speciesDB_compatHash(dbM, 298.15) != hB, "N2 nasa9_low[2] +0.001 -> different compat hash");
    check(speciesDB_compatHash(dbB, 0.0) != hB, "thermoHrefTemp 298.15 vs 0 -> different compat hash");
    check(speciesDB_compatHash(dbR, 298.15) != hB, "species order -> different compat hash");
    check(speciesDB_compatText(dbB, 298.15).find("builtin") == std::string::npos, "compat text has no source");

    // ---- (3) 記録 ----
    const fs::path runA = tmp / "runA", runM = tmp / "runM", runX = tmp / "runX";
    for (const auto& d : {runA, runM, runX}) fs::create_directories(d);
    const SpeciesRecordInfo rA = speciesDB_writeRecord(dbB, 298.15, "", "nozzle.h5", "verified", 0, runA.string());
    check(rA.recordFile == "resolved_species_" + hB.substr(0, 16) + ".yaml", "record name = resolved_species_<compat16>.yaml");
    check(speciesDB_sha256Hex(readAll(runA / rA.recordFile)) == rA.recordSha256, "record integrity = sha256(file)");
    const SpeciesRecordInfo rA2 = speciesDB_writeRecord(dbB, 298.15, "", "nozzle.h5", "verified", 0, runA.string());
    check(rA2.recordFile == rA.recordFile && rA2.recordSha256 == rA.recordSha256, "same record reused (not overwritten)");
    const SpeciesRecordInfo rA3 = speciesDB_writeRecord(dbB, 298.15, "", "nozzle.h5", "unverified_env", 1, runA.string());
    check(rA3.recordFile != rA.recordFile && rA3.compatHash == rA.compatHash, "different provenance -> separate record file, same compat hash: " + rA3.recordFile);
    check(readAll(runA / rA.recordFile).size() > 0 && speciesDB_sha256Hex(readAll(runA / rA.recordFile)) == rA.recordSha256,
          "original record untouched");
    {
        const auto d = speciesDB_diffRecord((runA / rA.recordFile).string(), dbB, 298.15);
        check(d.empty(), "diffRecord(record, same db) empty");
        const auto dm = speciesDB_diffRecord((runA / rA.recordFile).string(), dbM, 298.15);
        bool hit = false;
        for (const auto& x : dm) if (has(x, "N2.nasa9_low[2]")) hit = true;
        check(dm.size() == 1 && hit, "diffRecord shows only N2.nasa9_low[2]" + (dm.empty() ? std::string() : ": " + dm.front()));
    }
    const SpeciesRecordInfo rM = speciesDB_writeRecord(dbM, 298.15, dbMod.string(), "", "not_checked_resolve_only", 0, runM.string());

    // ---- (4) 入力場の照合 ----
    std::string st, msg; int unv = -1;
    const std::vector<std::string> dirsA = {runA.string()};
    check(speciesDB_checkInputField(dbB, 298.15, hB, rA.recordSha256, 0, "A.h5", dirsA, false, st, unv, msg)
          && st == "verified" && unv == 0, "A: same species -> allowed (verified)");
    check(speciesDB_checkInputField(dbF, 298.15, hB, rA.recordSha256, 0, "A.h5", dirsA, false, st, unv, msg),
          "source-only difference (external DB with identical coefficients) -> allowed");
    check(speciesDB_checkInputField(dbB, 298.15, hB, rA.recordSha256, 1, "A.h5", dirsA, false, st, unv, msg)
          && st == "unverified_inherited" && unv == 1, "unverified mark on input is inherited");
    // 属性なし
    check(!speciesDB_checkInputField(dbB, 298.15, "", "", -1, "old.h5", dirsA, false, st, unv, msg)
          && has(msg, "UNVERIFIABLE") && has(msg, "FORGE_ALLOW_UNVERIFIED_SPECIES=1"), "no attribute -> UNVERIFIABLE, env guidance");
    check(speciesDB_checkInputField(dbB, 298.15, "", "", -1, "old.h5", dirsA, true, st, unv, msg)
          && st == "unverified_env" && unv == 1, "no attribute + env -> allowed, marked unverified");
    // B/(c): 外部 DB で N2 low[2] +0.001 の run が A の場を読む → 拒否、該当係数を表示
    const bool okB = speciesDB_checkInputField(dbM, 298.15, hB, rA.recordSha256, 0, (runA / "res_100.h5").string(),
                                               {runA.string(), runM.string()}, false, st, unv, msg);
    check(!okB && has(msg, "N2.nasa9_low[2]"), "B/(c): modified N2 low[2] -> refused with coefficient shown");
    std::printf("---- message (B) ----\n%s\n---------------------\n", msg.c_str());
    check(!speciesDB_checkInputField(dbM, 298.15, hB, rA.recordSha256, 0, "x.h5", dirsA, true, st, unv, msg),
          "mismatch is not allowed by env");
    // 記録なし
    check(!speciesDB_checkInputField(dbM, 298.15, hB, rA.recordSha256, 0, "x.h5", {runX.string()}, false, st, unv, msg)
          && has(msg, "cannot be identified"), "mismatch without record -> 'cannot be identified'");
    // (b) 記録の取り違え: runX に runM の記録を A の名前で置く
    fs::copy_file(runM / rM.recordFile, runX / rA.recordFile);
    check(!speciesDB_checkInputField(dbM, 298.15, hB, rA.recordSha256, 0, "x.h5", {runX.string()}, false, st, unv, msg)
          && has(msg, "integrity"), "(b): swapped record -> integrity mismatch reported");
    {
        const auto d = speciesDB_diffRecord((runX / rA.recordFile).string(), dbB, 298.15);
        bool self = false;
        for (const auto& x : d) if (has(x, "compat_hash") || has(x, "nasa9_low[2]")) self = true;
        check(self, "(b): swapped record content differs from the name/species");
    }
    // 記録内容の改竄 (係数を書き換えて compat_hash はそのまま) → 自己整合の不一致
    {
        std::string t = readAll(runA / rA.recordFile);
        const std::string key = "nasa9_low: [";
        const size_t p = t.find(key);
        t.insert(p + key.size(), "1");   // 先頭係数の頭に 1 を足す
        std::ofstream(runX / "edited.yaml", std::ios::binary) << t;
        const auto d = speciesDB_diffRecord((runX / "edited.yaml").string(), dbB, 298.15);
        bool self = false;
        for (const auto& x : d) if (has(x, "recomputed")) self = true;
        check(self, "edited record -> compat_hash self-consistency mismatch");
    }

    fs::remove_all(tmp);
    std::printf("%s (%d failure%s)\n", g_fail ? "FAILED" : "ALL PASSED", g_fail, g_fail == 1 ? "" : "s");
    return g_fail ? 1 : 0;
}
