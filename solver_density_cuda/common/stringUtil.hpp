#ifndef STRING_UTIL_H
#define STRING_UTIL_H

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <cctype>
#include <algorithm>

// 空白文字 (space, \t, \n, \v, \f, \r) で分割する。
// 従来使っていた boost::split(out, in, boost::is_space()) と同じ規約:
//   - 連続する区切りは圧縮しない (token_compress_off) → 空トークンをそのまま残す
//   - 入力が空でも要素数 1 (空文字列) を返す
//   - 先頭・末尾の区切りも空トークンを生む (CRLF の \r は区切り扱い)
// gmsh の .msh は列位置でインデックスアクセスするため、この規約を崩すと静かに壊れる。
inline void splitOnSpace(std::vector<std::string> &out, const std::string &in)
{
    out.clear();
    std::string::size_type start = 0;
    for (std::string::size_type i = 0; i <= in.size(); ++i) {
        if (i == in.size() || std::isspace(static_cast<unsigned char>(in[i]))) {
            out.emplace_back(in, start, i - start);
            start = i + 1;
        }
    }
};

void eraseLastFeed(std::string &s)
{
    if (!s.empty() && s[s.length()-1] == '\n' && s[s.length()-1] == '\r') {
        s.erase(s.length()-1);
    };
};

#endif 