# =============================================================================
# embed_species_data.cmake — 共通 species データ (YAML テキスト) を C++ ヘッダへ埋め込む (cmake -P スクリプト)。
#   plans/active/thermophysics-solver-owned-species-db.md §5.1 #4。実行時のパス依存を避けるため、
#   ソルバ・変換器はこのヘッダ (生の文字列リテラル) を起動時に yaml-cpp でパースする (input/speciesDB.cpp)。
#
#   cmake -DIN=<forge_species_v1.yaml> -DOUT=<dir>/forge_species_data.hpp -P embed_species_data.cmake
#
# 出力:
#   kForgeSpeciesDataYaml[]    ファイル全文
#   kForgeSpeciesDataSha256[]  ファイル全文の SHA-256 (来歴用)
#   kForgeSpeciesDataName[]    ファイル名
# =============================================================================
if(NOT DEFINED IN OR NOT DEFINED OUT)
    message(FATAL_ERROR "embed_species_data.cmake: -DIN=<yaml> -DOUT=<header> are required")
endif()
if(NOT EXISTS "${IN}")
    message(FATAL_ERROR "embed_species_data.cmake: species data file not found: ${IN}")
endif()

file(READ "${IN}" _text)
file(SHA256 "${IN}" _sha)
get_filename_component(_name "${IN}" NAME)

set(_delim "FSPDATA")   # 生の文字列リテラルの区切りは 16 文字以内
string(FIND "${_text}" ")${_delim}\"" _pos)
if(NOT _pos EQUAL -1)
    message(FATAL_ERROR "embed_species_data.cmake: ${IN} contains the raw-string delimiter ')${_delim}\"'")
endif()

set(_hdr "// 生成ファイル (cmake/embed_species_data.cmake)。編集しない。元: ${_name}\n")
string(APPEND _hdr "#pragma once\n")
string(APPEND _hdr "static const char kForgeSpeciesDataName[] = \"${_name}\";\n")
string(APPEND _hdr "static const char kForgeSpeciesDataSha256[] = \"${_sha}\";\n")
string(APPEND _hdr "static const char kForgeSpeciesDataYaml[] = R\"${_delim}(${_text})${_delim}\";\n")

file(WRITE "${OUT}" "${_hdr}")
