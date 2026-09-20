#pragma once

#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <string>

#include "flowFormat.hpp"
#include "input/solverConfig.hpp"
#include "mesh/mesh.hpp"
#include "variables.hpp"

#include "yaml-cpp/yaml.h"

#include "cuda_forge/fluct_variables_d.cuh"


using namespace std;

struct bcondConfFormat{
    int physID;
    std::string physName;
    std::string kind;
    int outputHDFflg;
    std::map<std::string, int> inputInts ;
    std::map<std::string, flow_float> inputFloats;
    std::map<std::string, int> neuDirFlag ;

    map<string, map<string,int>> valueTypesOfBC
    {
        { // 0: float variables, 1: uniform float variables (read), 3: int variables , 
          //10:float information, 11:int information
          {"wall" , {
              {"ro"  ,0},
              {"roUx",0},
              {"roUy",0},
              {"roUz",0},
              {"roe" ,0},
              {"Ux"  ,1},
              {"Uy"  ,1},
              {"Uz"  ,1},
              {"Tt"  ,0},
              {"Pt"  ,0},
              {"Ts"  ,0},
              {"Ps"  ,0},
              {"kb"  ,0},
              {"omegab",0},
              {"ypls",0},
              {"twall_x",0},
              {"twall_y",0},
              {"twall_z",0},
              {"utau",0},
              {"qwall",0},   // WMLES 壁モデルの q_w (断熱では 0。viscousFlux wallTreatment==2 が消費)
          }},

          {"wall_isothermal", {
              {"ro"  ,0},
              {"roUx",0},
              {"roUy",0},
              {"roUz",0},
              {"roe" ,0},
              {"Ux"  ,1},
              {"Uy"  ,1},
              {"Uz"  ,1},
              {"Tt"  ,0},
              {"Pt"  ,0},
              {"Ts"  ,1},
              {"Ps"  ,0},
              {"kb"  ,0},
              {"omegab",0},
              {"ypls",0},
              {"twall_x",0},
              {"twall_y",0},
              {"twall_z",0},
              {"utau",0},
              {"qwall",0},   // WMLES 壁モデルの q_w (Kader。viscousFlux wallTreatment==2 が消費)
              // CHT の**保存的な実効界面熱量** $Q_f=\sum F^E-C$ の素材 (plan boundary-conjugate-heat-transfer §4.3)。
              //   `ifaceFw`   : 壁半割面が res_roe に入れた寄与そのもの (= $-\sum F^E_{\partial w}$)
              //   `ifaceRraw` : 壁ノードの**残差射影より前**の res_roe (= $R^{raw}$、定常 Dirichlet では $C=-R^{raw}$)
              // どちらも `output.interfaceDiag: 1` のときだけ書かれる (既定は触らないのでビット不変)。
              {"ifaceFw",0},
              {"ifaceRraw",0},

          }},

          {"inlet_uniformVelocity", {
              {"ro"  ,1},
              {"roUx",0},
              {"roUy",0},
              {"roUz",0},
              {"roe" ,0},
              {"Ux"  ,1},
              {"Uy"  ,1},
              {"Uz"  ,1},
              {"Tt"  ,0},
              {"Pt"  ,0},
              {"Ts"  ,0},
              {"Ps"  ,1},
              {"k"   ,1},
              {"omega",1},

          }},

          {"inlet_fluctVelocity", { 
              {"ro"  ,1},
              {"roUx",0},
              {"roUy",0},
              {"roUz",0},
              {"roe" ,0},
              {"Ux"  ,0},
              {"Uy"  ,0},
              {"Uz"  ,0},
              {"Ux0" ,1},
              {"Uy0" ,1},
              {"Uz0" ,1},
              {"Tt"  ,0},
              {"Pt"  ,0},
              {"Ts"  ,0},
              {"Ps"  ,1},
              {"k"   ,1},
              {"omega",1},

          }},


          {"inlet_Pressure", { 
              {"ro"  ,0},
              {"roUx",0},
              {"roUy",0},
              {"roUz",0},
              {"roe" ,0},
              {"Ux"  ,0},
              {"Uy"  ,0},
              {"Uz"  ,0},
              {"Tt"  ,1},
              {"Pt"  ,1},
              {"Ts"  ,0},
              {"Ps"  ,0},
              {"k"   ,1},
              {"omega",1},
          }},

          {"inlet_Pressure_dir", { 
              {"ro"  ,0},
              {"roUx",0},
              {"roUy",0},
              {"roUz",0},
              {"roe" ,0},
              {"Ux"  ,1},
              {"Uy"  ,1},
              {"Uz"  ,1},
              {"Tt"  ,1},
              {"Pt"  ,1},
              {"Ts"  ,0},
              {"Ps"  ,0},
              {"k"   ,1},
              {"omega",1},
          }},

          {"outlet_statPress", { 
              {"ro"  ,0},
              {"roUx",0},
              {"roUy",0},
              {"roUz",0},
              {"roe" ,0},
              {"Ux"  ,0},
              {"Uy"  ,0},
              {"Uz"  ,0},
              {"Tt"  ,1},
              {"Pt"  ,1},
              {"Ts"  ,0},
              {"Ps"  ,1},
              {"kb"  ,0},
              {"omegab",0},
          }},

          {"outflow", { 
              {"ro"  ,0},
              {"roUx",0},
              {"roUy",0},
              {"roUz",0},
              {"roe" ,0},
              {"Ux"  ,0},
              {"Uy"  ,0},
              {"Uz"  ,0},
              {"Tt"  ,1},
              {"Pt"  ,1},
              {"Ts"  ,0},
              {"Ps"  ,0},
              {"kb"  ,0},
              {"omegab",0},

          }},
 
          {"slip", { 
              {"ro"  ,-1},
              {"roUx",-1},
              {"roUy",-1},
              {"roUz",-1},
              {"roe" ,-1},
              {"Ux"  ,-1},
              {"Uy"  ,-1},
              {"Uz"  ,-1},
              {"Tt"  ,-1},
              {"Pt"  ,-1},
              {"Ts"  ,-1},
              {"Ps"  ,-1},
              {"kb"  ,-1},
              {"omegab",-1},
          }},

          {"axis", { 
              {"ro"  ,-1},
              {"roUx",-1},
              {"roUy",-1},
              {"roUz",-1},
              {"roe" ,-1},
              {"Ux"  ,-1},
              {"Uy"  ,-1},
              {"Uz"  ,-1},
              {"Tt"  ,-1},
              {"Pt"  ,-1},
              {"Ts"  ,-1},
              {"Ps"  ,-1},
              {"kb"  ,-1},
              {"omegab",-1},
          }},

          {"periodic", { 
              {"ro"  ,0},
              {"roUx",0},
              {"roUy",0},
              {"roUz",0},
              {"roe" ,0},
              {"Ux"  ,0},
              {"Uy"  ,0},
              {"Uz"  ,0},
              {"Tt"  ,0},
              {"Pt"  ,0},
              {"Ts"  ,0},
              {"Ps"  ,0},
              {"kb"  ,0},
              {"omegab",0},
              {"dtheta", 10},
              {"dx", 10},
              {"dy", 10},
              {"dz", 10},
              {"type" , 11},
              {"partnerBCID"  , 11},
              {"partnerPlnID" , 3},
              {"partnerCellID", 3},
          }},
        }
    };

    bcondConfFormat();
};

//void setBcondsValue(solverConfig& cfg , mesh& msh , variables& var , matrix& mat_p);
void applyBconds(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var , matrix& mat_p , fluct_variables& fluct);

void applyRansScalarBoundaries(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

void copyBcondsGradient(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var , matrix& mat_p);

void readBcondConfig(solverConfig& , vector<bcond>& );

// 入口分布プロファイル: inlet bcond の per-face bvar を CSV (`inlet_profile_<physID>.csv`) から
// face 重心座標で補間してセットする (ints:{inletProfile:1} で有効化)。readBcondConfig の後・
// 最初の applyBconds より前に呼ぶ。詳細は boundaryCond.cpp の関数ヘッダ参照。
void applyInletProfiles(solverConfig& cfg , mesh& msh);

// 壁温分布プロファイル: `ints: {wallProfile: 1}` の壁 bcond の per-face `Ts` を
// `wall_profile_<physID>.csv` から補間してセットする (inlet 版と同じ CSV 書式)。
// **評価点は値を課す位置** (node モードは壁ノード座標、cell モードは CV 重心) で、
// inlet の face 重心とは違う。readBcondConfig の後・最初の applyBconds より前に呼ぶ。
void applyWallProfiles(solverConfig& cfg , mesh& msh);

