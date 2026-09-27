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
#   kForgeTransportDataYaml[] / kForgeTransportDataSha256[] / kForgeTransportDataName[]
#                              同じディレクトリの輸送データ forge_transport_v1.yaml (CEA trans.inp から生成;
#                              tools/cea_trans_to_forge_transport.py、plan #5t2)。IN と同じ場所に必須 (-DTRANS で上書き可)。
#                              同じヘッダに入れるのは、既存の単体試験の生成手順 (-DIN だけを渡す) をそのまま使うため。
# =============================================================================
if(NOT DEFINED IN OR NOT DEFINED OUT)
    message(FATAL_ERROR "embed_species_data.cmake: -DIN=<yaml> -DOUT=<header> are required")
endif()
if(NOT EXISTS "${IN}")
    message(FATAL_ERROR "embed_species_data.cmake: species data file not found: ${IN}")
endif()

if(NOT DEFINED TRANS)
    get_filename_component(_dir "${IN}" DIRECTORY)
    set(TRANS "${_dir}/forge_transport_v1.yaml")
endif()
if(NOT EXISTS "${TRANS}")
    message(FATAL_ERROR "embed_species_data.cmake: transport data file not found: ${TRANS}")
endif()

file(READ "${IN}" _text)
file(SHA256 "${IN}" _sha)
get_filename_component(_name "${IN}" NAME)

file(READ "${TRANS}" _ttext)
file(SHA256 "${TRANS}" _tsha)
get_filename_component(_tname "${TRANS}" NAME)

set(_delim "FSPDATA")   # 生の文字列リテラルの区切りは 16 文字以内
foreach(_pair "IN;_text" "TRANS;_ttext")
    list(GET _pair 0 _f)
    list(GET _pair 1 _v)
    string(FIND "${${_v}}" ")${_delim}\"" _pos)
    if(NOT _pos EQUAL -1)
        message(FATAL_ERROR "embed_species_data.cmake: ${${_f}} contains the raw-string delimiter ')${_delim}\"'")
    endif()
endforeach()

set(_hdr "// 生成ファイル (cmake/embed_species_data.cmake)。編集しない。元: ${_name}\n")
string(APPEND _hdr "#pragma once\n")
string(APPEND _hdr "static const char kForgeSpeciesDataName[] = \"${_name}\";\n")
string(APPEND _hdr "static const char kForgeSpeciesDataSha256[] = \"${_sha}\";\n")
string(APPEND _hdr "static const char kForgeSpeciesDataYaml[] = R\"${_delim}(${_text})${_delim}\";\n")
string(APPEND _hdr "static const char kForgeTransportDataName[] = \"${_tname}\";\n")
string(APPEND _hdr "static const char kForgeTransportDataSha256[] = \"${_tsha}\";\n")
string(APPEND _hdr "static const char kForgeTransportDataYaml[] = R\"${_delim}(${_ttext})${_delim}\";\n")

file(WRITE "${OUT}" "${_hdr}")
