#pragma once

#include "mesh/mesh.hpp"
#include "variables.hpp"
#include "solverConfig.hpp"
#include "calcWallDistance_kdtree.hpp"

#include <cstdlib>
#include <cmath>

void setInitial(solverConfig& cfg , mesh& msh , variables& v)
{
    std::cout << "cfg initial=" << cfg.initial << "\n";

    // *** set initial value ***
    // 1D_shock_tube
    if (cfg.initial == "sod") {
        flow_float gam  = 1.4;
        flow_float roL  = 1.0;
        flow_float prsL = 1.0;
        flow_float velL = 0.0;
        flow_float roR  = 0.125;
        flow_float prsR = 0.1;
        flow_float velR = 0.0;

        for (geom_int i = 0 ; i<msh.nCells ; i++) {
            if (msh.cells[i].centCoords[0] < 0.5) {
                v.c["ro"][i] = roL;
                v.c["roUx"][i] = roL*velL;
                v.c["roUy"][i] = 0.0;
                v.c["roUz"][i] = 0.0;
                v.c["roe"][i] = prsL/(gam-1.0) + 0.5*roL*velL*velL;
            } else {
                v.c["ro"][i] = roR;
                v.c["roUx"][i] = roR*velR;
                v.c["roUy"][i] = 0.0;
                v.c["roUz"][i] = 0.0;
                v.c["roe"][i] = prsR/(gam-1.0) + 0.5*roR*velR*velR;
            }
        }
    } else if (cfg.initial == "mach3") {
        // mach3 
        flow_float gam  = 1.4;
        flow_float roL  = 1.4;
        flow_float prsL = 1.0;
        flow_float velL = 3.0;

        for (geom_int i = 0 ; i<msh.nCells ; i++) {
            v.c["ro"][i] = roL;
            v.c["roUx"][i] = roL*velL;
            v.c["roUy"][i] = 0.0;
            v.c["roUz"][i] = 0.0;
            v.c["roe"][i] = prsL/(gam-1.0) + 0.5*roL*velL*velL;
        }

    } else if (cfg.initial == "bump") {
        // bump
        flow_float gam  = 1.4;
        flow_float Pt = 100000.0;
        flow_float Tt = 293.15 ;
        flow_float M  = 0.675;
        //flow_float M  = 1.65;

        flow_float Ps = Pt*pow((1.0+0.5*(gam-1.0)*M*M), -gam/(gam-1.0));
        flow_float Ts = Tt*pow((1.0+0.5*(gam-1.0)*M*M), -1.0);
        flow_float ro = gam*Ps/(cfg.cp*(gam-1.0)*Ts);
        flow_float a = sqrt(gam*Ps/ro);

        for (geom_int i = 0 ; i<msh.nCells ; i++) {

            v.c["ro"][i]   = ro;
            v.c["roUx"][i] = M*a*ro;
            v.c["roUy"][i] = 0.0*ro;
            v.c["roUz"][i] = 0.0*ro;
            v.c["roe"][i] =  pow(1.0+0.5*(gam-1.0)*M*M,-gam/(gam-1.0))*Pt/(gam-1.0) + 0.5*ro*pow((M*a),2.0);
        }


    } else if (cfg.initial == "bump_M1.65") {
        // bump
        flow_float gam  = 1.4;
        flow_float Pt = 100000.0;
        flow_float Tt = 293.15 ;
        //flow_float M  = 0.675;
        flow_float M  = 1.65;

        flow_float Ps = Pt*pow((1.0+0.5*(gam-1.0)*M*M), -gam/(gam-1.0));
        flow_float Ts = Tt*pow((1.0+0.5*(gam-1.0)*M*M), -1.0);
        flow_float ro = gam*Ps/(cfg.cp*(gam-1.0)*Ts);
        flow_float a = sqrt(gam*Ps/ro);


        for (geom_int i = 0 ; i<msh.nCells ; i++) {

            v.c["ro"][i]   = ro;
            v.c["roUx"][i] = M*a*ro;
            v.c["roUy"][i] = 0.0*ro;
            v.c["roUz"][i] = 0.0*ro;
            v.c["roe"][i] =  pow(1.0+0.5*(gam-1.0)*M*M,-gam/(gam-1.0))*Pt/(gam-1.0) + 0.5*ro*pow((M*a),2.0);
        }

    } else if (cfg.initial == "Taylor-Green" or cfg.initial == "Taylor-Green_M0.1") {
        flow_float gam = 1.4;

        flow_float L   = 1.0; // length

        flow_float ro0 = 1.0;
        flow_float P0  = 1.0/gam;
        flow_float M0  = 0.4;

        if (cfg.initial == "Taylor-Green_M0.1") M0 = 0.1;

        for (geom_int i = 0 ; i<msh.nCells ; i++) {

            flow_float x = msh.cells[i].centCoords[0];
            flow_float y = msh.cells[i].centCoords[1];
            flow_float z = msh.cells[i].centCoords[2];

            flow_float u0 = M0*sin(x/L)*cos(y/L)*cos(z/L);
            flow_float v0 =-M0*cos(x/L)*sin(y/L)*cos(z/L);
            flow_float w0 = 0.0;

            flow_float P1 = P0 +ro0*M0*M0/16.0*(cos(2*x/L)+cos(2*y/L))*(cos(2*z/L)+2.0);

            //flow_float ro1 = P1/(R*T0);

            v.c["ro"][i]   = ro0;
            v.c["roUx"][i] = ro0*u0;
            v.c["roUy"][i] = ro0*v0;
            v.c["roUz"][i] = 0.0;
            v.c["roe"][i]  = P1/(gam-1.0) + 0.5*ro0*(u0*u0+v0*v0+w0*w0);
        }

    } else if (cfg.initial == "half_sphere") {
        flow_float cp  = cfg.cp;
        flow_float gam = cfg.gamma;
        flow_float R   = cp*(gam-1.0)/gam;

        flow_float ro0 = 0.02026;
        flow_float T0  = 63.73;
        flow_float P0  = ro0*R*T0;

        flow_float u0  = 1296.22;

        for (geom_int i = 0 ; i<msh.nCells ; i++) {
            v.c["ro"][i]   = ro0;
            v.c["roUx"][i] = ro0*u0;
            v.c["roUy"][i] = 0.0;
            v.c["roUz"][i] = 0.0;
            v.c["roe"][i]  = P0/(gam-1.0) + 0.5*ro0*(u0*u0);
        }

    } else if (cfg.initial == "supercritical") {
        flow_float cp  = cfg.cp;
        flow_float gam = cfg.gamma;
        flow_float R   = cp*(gam-1.0)/gam;

        flow_float M = 0.3;
        flow_float P = 101325.0;
        flow_float T = 288.15;

        flow_float ro= P/(R*T);
        flow_float c = sqrt(gam*R*T);
        flow_float u = c*M;

        for (geom_int i = 0 ; i<msh.nCells ; i++) {

            v.c["ro"][i]   = ro;
            v.c["roUx"][i] = ro*u;
            v.c["roUy"][i] = 0.0;
            v.c["roUz"][i] = 0.0;
            v.c["roe"][i]  = P/(gam-1.0) + 0.5*ro*(u*u);
        }

    } else if (cfg.initial == "uniform_p101325_u10") {
        flow_float cp  = cfg.cp;
        flow_float gam = cfg.gamma;
        flow_float R   = cp*(gam-1.0)/gam;

        flow_float P = 101325.0;
        flow_float T = 300.0;
        flow_float u = 10.0;

        flow_float ro = P/(R*T);

        for (geom_int i = 0 ; i<msh.nCells ; i++) {

            v.c["ro"][i]   = ro;
            v.c["roUx"][i] = ro*u;
            v.c["roUy"][i] = 0.0;
            v.c["roUz"][i] = 0.0;
            v.c["roe"][i]  = P/(gam-1.0) + 0.5*ro*(u*u);
        }

    } else if (cfg.initial == "slot_m2") {
        // case/58 深いスロットの CHT (V6′, plan boundary-conjugate-heat-transfer §6)。
        // 入口 inlet_uniformVelocity (ro=0.116144, U=694.38, Ps=10000) と**同じ状態**を一様に置いて
        // 起動過渡を抑える。**`uniform_p101325_u10` を使ってはいけない**: P が 10 倍・U が 1/70 で入口と
        // 不整合なため、2000 step で P が 2.18e7 Pa (自由流の 2000 倍)・ro 175 kg/m3 に達し
        // `rms_roUy` が RISING になった (2026-09-26 実測。AGENTS.md「発散の主因 = 初期値が入口流れと不整合」)。
        // 空気 CPG (gamma=1.4), Ts=300 K, M=2.0。袋小路のスロットも同じ状態で満たす。
        flow_float cp  = cfg.cp;
        flow_float gam = cfg.gamma;
        flow_float R   = cp*(gam-1.0)/gam;

        flow_float Ps = 10000.0;
        flow_float Ts = 300.0;
        flow_float M  = 2.0;

        flow_float ro = Ps/(R*Ts);
        flow_float c  = sqrt(gam*R*Ts);
        flow_float u  = c*M;

        for (geom_int i = 0 ; i<msh.nCells ; i++) {

            v.c["ro"][i]   = ro;
            v.c["roUx"][i] = ro*u;
            v.c["roUy"][i] = 0.0;
            v.c["roUz"][i] = 0.0;
            v.c["roe"][i]  = Ps/(gam-1.0) + 0.5*ro*(u*u);
        }

    } else if (cfg.initial == "laval") {
        flow_float cp  = cfg.cp;
        flow_float gam = cfg.gamma;
        flow_float R   = cp*(gam-1.0)/gam;

        flow_float M = 0.80;
        flow_float P0 = 59100.0;
        flow_float T0 = 286.65;

        flow_float Ps = P0/(pow(1.0+0.5*(gam-1.0)*M*M, gam/(gam-1.0)));
        flow_float Ts = T0/(1.0+0.5*(gam-1.0)*M*M);

        flow_float ro= Ps/(R*Ts);
        flow_float c = sqrt(gam*R*Ts);
        flow_float u = c*M;

        for (geom_int i = 0 ; i<msh.nCells ; i++) {
            v.c["ro"][i]   = ro;
            v.c["roUx"][i] = ro*u;
            v.c["roUy"][i] = 0.0;
            v.c["roUz"][i] = 0.0;
            v.c["roe"][i]  = Ps/(gam-1.0) + 0.5*ro*(u*u);
        }

    } else if (cfg.initial == "hifire") {
        flow_float cp  = cfg.cp;
        flow_float gam = cfg.gamma;
        flow_float R   = cp*(gam-1.0)/gam;

        flow_float M = 7.16;
        flow_float P = 4620.0;
        flow_float T = 231.7;

        flow_float ro= P/(R*T);
        flow_float c = sqrt(gam*R*T);
        flow_float u = c*M;

        for (geom_int i = 0 ; i<msh.nCells ; i++) {
            v.c["ro"][i]   = ro;
            v.c["roUx"][i] = ro*u;
            v.c["roUy"][i] = 0.0;
            v.c["roUz"][i] = 0.0;
            v.c["roe"][i]  = P/(gam-1.0) + 0.5*ro*(u*u);
        }

    } else if (cfg.initial == "flare") {
        flow_float cp  = cfg.cp;
        flow_float gam = cfg.gamma;
        flow_float R   = cp*(gam-1.0)/gam;

        flow_float M = 6.0;
        //flow_float P = 930.9*1000;
        //flow_float T = 433.0;
        flow_float T = 52.805;

        flow_float ro= 0.03888;
        flow_float P= ro*R*T;
        flow_float c = sqrt(gam*R*T);
        flow_float u = c*M;

        for (geom_int i = 0 ; i<msh.nCells ; i++) {
            v.c["ro"][i]   = ro;
            v.c["roUx"][i] = ro*u;
            v.c["roUy"][i] = 0.0;
            v.c["roUz"][i] = 0.0;
            v.c["roe"][i]  = P/(gam-1.0) + 0.5*ro*(u*u);
        }

    } else if (cfg.initial == "poiseuille") {
        flow_float cp  = cfg.cp;
        flow_float gam = cfg.gamma;
        flow_float R   = cp*(gam-1.0)/gam;

        flow_float P = 100000.0;
        flow_float T = 288.15;

        flow_float ro= P/(R*T);
        flow_float u = 10;

        for (geom_int i = 0 ; i<msh.nCells ; i++) {
            v.c["ro"][i]   = ro;
            v.c["roUx"][i] = ro*u;
            v.c["roUy"][i] = 0.0;
            v.c["roUz"][i] = 0.0;
            v.c["roe"][i]  = P/(gam-1.0) + 0.5*ro*(u*u);
        }

    } else if (cfg.initial == "backstep") {
        // bump
        flow_float gam  = 1.4;
        flow_float Pt = 100000.0;
        flow_float Tt = 288.15 ;
        flow_float M  = 0.1 ;

        flow_float Ps = Pt*pow((1.0+0.5*(gam-1.0)*M*M), -gam/(gam-1.0));
        flow_float Ts = Tt*pow((1.0+0.5*(gam-1.0)*M*M), -1.0);
        flow_float ro = gam*Ps/(cfg.cp*(gam-1.0)*Ts);
        flow_float a = sqrt(gam*Ps/ro);

        flow_float u = M*a;

        // RANS (SST) 用の自由流乱流量シード。冷間始動で roK=roOmega=0 だと
        // mut=ro*k/omega=0 のまま乱流が立ち上がらないため、freestream 値を与える。
        // mu_t/mu ~ 10 を目安に k=4.0, omega=500.0 (inlet bcond の k/omega と整合)。
        // LES (WALE) では roK/roOmega は未使用なので無害。
        flow_float kIni     = 4.0;
        flow_float omegaIni = 500.0;

        for (geom_int i = 0 ; i<msh.nCells ; i++) {

            v.c["ro"][i]   = ro;
            v.c["roUx"][i] = M*a*ro;
            v.c["roUy"][i] = 0.0*ro;
            v.c["roUz"][i] = 0.0*ro;
            //v.c["roe"][i] =  pow(1.0+0.5*(gam-1.0)*M*M,-gam/(gam-1.0))*Pt/(gam-1.0) + 0.5*ro*pow((M*a),2.0);
            v.c["roe"][i]  = Ps/(gam-1.0) + 0.5*ro*(u*u);
            v.c["roK"][i]     = ro*kIni;
            v.c["roOmega"][i] = ro*omegaIni;
        }

    } else if (cfg.initial == "flat_plate") {
        // 乱流平板 (SST 壁法則検証)。入口全圧/全温と M から自由流の静的状態を決め、
        // 全域を一様流で初期化する。粘性は realistic air (mu~1.8e-5)。
        // freestream 乱流量が小さすぎると壁近傍で k->0 に固定され SST が相対層流化する
        // (Pk = mu_t S^2, mu_t=rho a1 k/omega なので k=0 は不動点)。
        // そこで実績のある naca_ml RANS と同等の乱流レベル (mu_t/mu ~ 60) を与え、
        // かつ omega を下げて長い領域での自由流減衰を抑える (k=0.3, omega=300)。
        // inlet bcond の k/omega と一致させること。
        flow_float gam  = 1.4;
        flow_float Pt = 100000.0;
        flow_float Tt = 288.15 ;
        flow_float M  = 0.2 ;

        flow_float Ps = Pt*pow((1.0+0.5*(gam-1.0)*M*M), -gam/(gam-1.0));
        flow_float Ts = Tt*pow((1.0+0.5*(gam-1.0)*M*M), -1.0);
        flow_float ro = gam*Ps/(cfg.cp*(gam-1.0)*Ts);
        flow_float a = sqrt(gam*Ps/ro);

        flow_float u = M*a;

        // freestream 乱流量: mu_t/mu ~ 60 (naca_ml 相当), omega は減衰抑制のため低め
        flow_float kIni     = 0.3;
        flow_float omegaIni = 300.0;

        for (geom_int i = 0 ; i<msh.nCells ; i++) {
            v.c["ro"][i]   = ro;
            v.c["roUx"][i] = u*ro;
            v.c["roUy"][i] = 0.0;
            v.c["roUz"][i] = 0.0;
            v.c["roe"][i]  = Ps/(gam-1.0) + 0.5*ro*(u*u);
            v.c["roK"][i]     = ro*kIni;
            v.c["roOmega"][i] = ro*omegaIni;
        }

    } else if (cfg.initial == "kusabi") {
        flow_float cp  = cfg.cp;
        flow_float gam = cfg.gamma;
        flow_float R   = cp*(gam-1.0)/gam;

        flow_float M = 6.0;
        flow_float P = 2500.0;
        flow_float T = 220.0;

        flow_float ro= P/(R*T);
        flow_float c = sqrt(gam*R*T);
        flow_float u = c*M;

        for (geom_int i = 0 ; i<msh.nCells ; i++) {
            v.c["ro"][i]   = ro;
            v.c["roUx"][i] = ro*u;
            v.c["roUy"][i] = 0.0;
            v.c["roUz"][i] = 0.0;
            v.c["roe"][i]  = P/(gam-1.0) + 0.5*ro*(u*u);
        }

    } else if (cfg.initial == "nozzle_wys") {
        // Wyslouzil supersonic nozzle 用の一様亜音速初期場。
        // 入口 inlet_Pressure (Pt=101325, Tt=293.15) と総圧をそろえ、起動過渡を抑える。
        flow_float cp  = cfg.cp;
        flow_float gam = cfg.gamma;
        flow_float R   = cp*(gam-1.0)/gam;

        flow_float M  = 0.30;
        flow_float P0 = 101325.0;
        flow_float T0 = 293.15;

        flow_float Ps = P0/(pow(1.0+0.5*(gam-1.0)*M*M, gam/(gam-1.0)));
        flow_float Ts = T0/(1.0+0.5*(gam-1.0)*M*M);

        flow_float ro = Ps/(R*Ts);
        flow_float c  = sqrt(gam*R*Ts);
        flow_float u  = c*M;

        for (geom_int i = 0 ; i<msh.nCells ; i++) {
            v.c["ro"][i]   = ro;
            v.c["roUx"][i] = ro*u;
            v.c["roUy"][i] = 0.0;
            v.c["roUz"][i] = 0.0;
            v.c["roe"][i]  = Ps/(gam-1.0) + 0.5*ro*(u*u);
        }

    } else if (cfg.initial == "arthur_n2") {
        // Arthur (1952) 2D source-flow ノズルの dry 超音速膨張用一様初期場。
        // 発散部のみを切り出し、入口 (スロート) の超音速状態 M=1.05 を一様に置く。
        // 入口 inlet_uniformVelocity (ρ,U,Ps 固定) と同じ状態にして起動過渡を抑える。
        // 貯気: P0=844037 Pa, T0=290 K (Arthur Run 34-1 / Fig.4,11 の P0=8.33 atm 相当)。N2: γ=1.4。
        flow_float gam = cfg.gamma;
        flow_float R   = cfg.cp*(gam-1.0)/gam;

        flow_float M  = 1.05;
        flow_float P0 = 844037.0;
        flow_float T0 = 290.0;

        flow_float Ts = T0/(1.0+0.5*(gam-1.0)*M*M);
        flow_float Ps = P0*pow(Ts/T0, gam/(gam-1.0));
        flow_float ro = Ps/(R*Ts);
        flow_float c  = sqrt(gam*R*Ts);
        flow_float u  = c*M;

        for (geom_int i = 0 ; i<msh.nCells ; i++) {
            v.c["ro"][i]   = ro;
            v.c["roUx"][i] = ro*u;
            v.c["roUy"][i] = 0.0;
            v.c["roUz"][i] = 0.0;
            v.c["roe"][i]  = Ps/(gam-1.0) + 0.5*ro*(u*u);
        }

    } else if (cfg.initial == "passive_pseudoshock") {
        // case/36 多孔壁パッシブコントロール (Matsuo et al. 1988) の超音速吹き抜け起動用。
        // 入口 (助走領域) の超音速状態 M=1.689 を一様に置き、inlet_uniformVelocity
        // (ρ,U,Ps 固定) と同じ状態で起動過渡を抑える。
        // 貯気: P0=3 MPa, T0=288.15 K (空気 CPG γ=1.4)。面積比から M_in=1.689。
        flow_float gam = cfg.gamma;
        flow_float R   = cfg.cp*(gam-1.0)/gam;

        flow_float M  = 1.689;
        flow_float P0 = 3.0e6;
        flow_float T0 = 288.15;

        flow_float Ts = T0/(1.0+0.5*(gam-1.0)*M*M);
        flow_float Ps = P0*pow(Ts/T0, gam/(gam-1.0));
        flow_float ro = Ps/(R*Ts);
        flow_float c  = sqrt(gam*R*Ts);
        flow_float u  = c*M;

        for (geom_int i = 0 ; i<msh.nCells ; i++) {
            v.c["ro"][i]   = ro;
            v.c["roUx"][i] = ro*u;
            v.c["roUy"][i] = 0.0;
            v.c["roUz"][i] = 0.0;
            v.c["roe"][i]  = Ps/(gam-1.0) + 0.5*ro*(u*u);
        }

    } else {
        std::cout << "Error: Unknown initial" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    std::list<std::string> names = {"ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega"};
    v.copyVariables_cell_H2D(names);

    calcWallDistance_kdtree(cfg , msh , v);
};


