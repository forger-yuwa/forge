#pragma once
// 変換器 (convertGmshToForge) の工程別メモリ計測 (plan tooling-sern-mesh-blocking §5.1 B4-5 (1))。
//
// 環境変数 FORGE_MEMLOG=1 のときだけ、工程の境目で /proc/self/status の VmRSS (現在) と VmHWM (ここまでのピーク) を
// 1 行ずつ stdout に出す。既定 (未設定・0) では何も出さず、出力 h5 にも触れない (計測だけ)。
// 主要コンテナの推定バイトは size/capacity × sizeof と glibc malloc のチャンク (最小 32 B, 16 B 境界) で概算する
// (ネストした vector の 1 要素ごとのヘッダ・ヒープ確保を含めるのが目的。厳密値ではない)。
//
// 使い方: MEMLOG("label", 追加文字列式)。追加文字列式は有効時だけ評価される (無効時は走査コストも無い)。

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <list>
#include <sstream>
#include <string>
#include <vector>

namespace memlog {

inline bool enabled()
{
    static const bool on = [] {
        const char* e = std::getenv("FORGE_MEMLOG");
        return e != nullptr && std::atoi(e) != 0;
    }();
    return on;
}

// /proc/self/status の VmRSS / VmHWM [kB]。読めなければ -1。
inline void readStatus(long& rssKB, long& hwmKB)
{
    rssKB = -1; hwmKB = -1;
    FILE* f = std::fopen("/proc/self/status", "r");
    if (!f) return;
    char buf[256];
    while (std::fgets(buf, sizeof(buf), f)) {
        if (std::strncmp(buf, "VmRSS:", 6) == 0) rssKB = std::atol(buf + 6);
        else if (std::strncmp(buf, "VmHWM:", 6) == 0) hwmKB = std::atol(buf + 6);
    }
    std::fclose(f);
}

// glibc malloc の 1 確保あたりの実消費 (ヘッダ 8 B + 16 B 境界、最小 32 B) の概算。0 B 要求は確保なしとみなす。
inline size_t heapChunk(size_t req)
{
    if (req == 0) return 0;
    size_t c = (req + 8 + 15) & ~size_t(15);
    return c < 32 ? 32 : c;
}

// 平坦な vector<T> の推定バイト (capacity × sizeof)。
template <class T>
inline size_t flatBytes(const std::vector<T>& v) { return heapChunk(v.capacity() * sizeof(T)); }

// vector<内部に vector を持つ構造体> の推定バイト。inner(e) が要素 e の内側ヒープ推定を返す。
template <class T, class F>
inline size_t nestedBytes(const std::vector<T>& v, F inner)
{
    size_t b = heapChunk(v.capacity() * sizeof(T));
    for (const auto& e : v) b += inner(e);
    return b;
}

template <class T>
inline size_t innerVec(const std::vector<T>& v) { return heapChunk(v.capacity() * sizeof(T)); }

// "名前[要素数]=xxx MB" 形式の 1 項目。
inline std::string item(const std::string& name, size_t n, size_t bytes)
{
    std::ostringstream os;
    os << name << "[" << n << "]=" << (double)bytes / 1048576.0 << "MB";
    return os.str();
}

inline void log(const char* label, const char* func, int line, const std::string& extra)
{
    long rss, hwm;
    readStatus(rss, hwm);
    std::cout << "[memlog] " << label << " @" << func << ":" << line
              << " VmRSS=" << rss / 1024 << "MB VmHWM=" << hwm / 1024 << "MB";
    if (!extra.empty()) std::cout << " | " << extra;
    std::cout << std::endl;
}

} // namespace memlog

#define MEMLOG(label, extra) \
    do { if (memlog::enabled()) memlog::log((label), __func__, __LINE__, (extra)); } while (0)
