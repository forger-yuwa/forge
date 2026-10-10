// =============================================================================
// 対流フラックス KEEP_d の実装断片。
//   - 純粋 KEEP (Kinetic Energy and Entropy Preserving) 中心流束。**中心部は散逸項なし**。
//     隣接セル/ノード対 (ic0,ic1) の生値で中心流束を構成する。KE は二次分割で保存、エントロピーは
//     **厳密保存ではない** (算術/交差平均構成のため Tadmor 条件 Δwᵀ·F* = Δψ を満たさない。
//     誤差はジャンプ 3 次 O(Δ³) の準保存 = 滑らかな場では設計次数で整合、大ジャンプでは非保存。
//     tools/verify_keep_tadmor.py で数値確認済。厳密 EC には log-mean 平均 (KEPEC 型) が必要)。
//     したがって散逸レイヤ込みでも「流束全体の entropy stability」は主張しない (散逸項単体が ES)。
//     Roe 行列散逸・MUSCL 再構成・リミタ・Ducros は持たない (低散逸 LES/ILES 用; SGS は WALE が担う)。
//   - opt-in 散逸レイヤ (keepDissType, plans/accepted/convection-keep-es-dissipation.md):
//     keepDissType==0 で従来どおり散逸ゼロ (ビット不変)。==1 で scalar entropy-stable 散逸
//     F -= 0.5*σ*λ'*ΔU を全 5 式に加算 (λ'=|Un|+c'; lowMachPrecond>=1 なら c'=lowMachCprime 低マッハ
//     スケール、else c)。LLF 型は Δw·ΔU>=0 (エントロピー凸性) で ES。純 KEEP は低マッハ圧力 2Δ 市松を
//     原理的に減衰できない (中心平均の null-mode) ため、その減衰はこのレイヤが担う。
//   - convectiveFlux_d.cu に textual include される断片であり **スタンドアロン TU ではない**。
//     CMake の source には足さないこと。
//   - 周回面は geom.nLoopPlanes (= convPlaneBound)。cell では内部+境界 ghost を周回し、境界面は
//     ic1=ghost セルの生値で中心流束を組む。node 弱形式では内部双対面のみ (境界は別カーネル/gather)。
// =============================================================================

__global__ void KEEP_d
(
 flow_float ga,
 int thermalMethod, const SpeciesThermo* sp,      // Step 3/4: TP (thermalMethod==2)。CPG では未使用
 int nSpecies, flow_float** roY,                  // Step 4: 多成分の組成 (nSpecies<=1 なら未参照)
 int keepDissType, flow_float keepDissCoeff,      // opt-in ES 散逸レイヤ (0: off=ビット不変)
 int keepDissCprime, flow_float precondEps,       // 散逸波速: 1 で音響 c'=lowMachCprime (lowMachPrecond から独立)
 int keepDissJump,                                // 散逸ジャンプ: 0=生 (既定) / 1=再構成後 (matrix CPG 枝のみ)
 int keepDissPrecond,                             // 1: 音響対散逸を Turkel 前処理 2×2 に置換 (matrix CPG 枝のみ)
 const flow_float* fd_shield,                     // f_d 駆動 σ ブレンド (nullptr=off, turbulence-iddes-sst §4.8)
 int fdLesOne,                                    //   fd_shield の向き: 1=DDES (f_d, 1=LES) / 0=IDDES (f̃_d, 1=RANS)
 flow_float keepDissCoeffMax,                     //   σ_f=max(keepDissCoeff, ransFrac·keepDissCoeffMax)
 flow_float keepDissCbCoeff,                      // 高周波圧力欠陥駆動 mass-flux 補正 C_cb (0=off・ビット不変)
 flow_float keepDissCbEps,                        //   その基準速度カットオフ ε_cb (Ur=min(c,max(|u|,ε·c)))
 int keepDissOpBlendRaw,                          // 1: opBlend 標準 Roe 側に生ジャンプを使う (RANS/grey 帯=完全風上相当)
 GradFields    grd,
 FaceGeom      geom,
 PrimState     st,
 ResidualOut   reso
)
{
    // --- struct 引数をローカルへ展開 ---
    const geom_int nCells       = geom.nCells;
    geom_int*      plane_cells  = geom.plane_cells;
    const geom_int nLoopPlanes  = geom.nLoopPlanes;
    geom_int*      loop_planes  = geom.loop_planes;
    geom_float *sx=geom.sx, *sy=geom.sy, *sz=geom.sz, *ss=geom.ss;
    flow_float* massflux=geom.massflux;

    flow_float *ro=st.ro, *Ux=st.Ux, *Uy=st.Uy, *Uz=st.Uz, *Ps=st.Ps, *roe_c=st.roe;

    flow_float *res_ro=reso.res_ro, *res_roUx=reso.res_roUx, *res_roUy=reso.res_roUy, *res_roUz=reso.res_roUz, *res_roe=reso.res_roe;
    // --- ローカル展開ここまで ---

    geom_int ip_orig = blockDim.x*blockIdx.x + threadIdx.x;

    if (ip_orig < nLoopPlanes) {

        geom_int ip = loop_planes[ip_orig];

        geom_int ic0 = plane_cells[2*ip+0];
        geom_int ic1 = plane_cells[2*ip+1];
        (void)nCells;

        geom_float sxx = sx[ip], syy = sy[ip], szz = sz[ip], sss = ss[ip];
        geom_float nx = sxx/sss, ny = syy/sss, nz = szz/sss;

        // ---- 純粋 KEEP 中心流束 (隣接対の生値; KE 保存・エントロピー O(Δ³) 準保存・散逸なし) ----
        // U∞≠0 の動く一様流保存 (advGauge, plan §8.2): roRef>0 (CPG のみ) では中心流束を
        // 「元の流束 − 参照一様流束 F∞(s)」の差分形因数分解 ((差分)×(平均)+(参照)×(差分)) で組む。
        // F∞ は s に線形な定数流束: 内部面で等価逆符号のまま=保存厳密、セル和は F∞(Σs)=0 の解析ゲージ。
        // 一様場 (ρ,u,p)=(ρ∞,u∞,p∞) では全因子がビット単位ゼロ → free-stream を機械精度で保存。
        // 式は tools/verify_advective_gauge.py で (1)解析等価 (2)一様 float32 ビット零
        // (3)orig−gauged==F∞(r_c) の telescoping を数値検証済。roRef=0 (既定) は従来経路でビット不変。
        const bool advGauge = (d_roRef > 0.0) && (thermalMethod != 2);
        flow_float res_ro_temp, res_roe_temp, res_roUx_temp, res_roUy_temp, res_roUz_temp;
        flow_float CinfS = 0.0;   // 参照質量流束 F∞·s (massflux の物理化用)

        if (!advGauge) {
        flow_float Ctilde  = 0.5*(ro[ic0]+ro[ic1])*0.5*( (Ux[ic0]+Ux[ic1])*nx
                                                        +(Uy[ic0]+Uy[ic1])*ny
                                                        +(Uz[ic0]+Uz[ic1])*nz );
        flow_float Mtildex = Ctilde*(Ux[ic0]+Ux[ic1])*0.5;
        flow_float Mtildey = Ctilde*(Uy[ic0]+Uy[ic1])*0.5;
        flow_float Mtildez = Ctilde*(Uz[ic0]+Uz[ic1])*0.5;
        flow_float Ktilde  = Ctilde*0.5*(Ux[ic0]*Ux[ic1] +Uy[ic0]*Uy[ic1] +Uz[ic0]*Uz[ic1]);
        flow_float Itilde;
        if (thermalMethod == 2) {
            // TP (Step 3): 内部エネルギーを保存量 roe から取る (EOS 非依存・保存整合)。
            // CPG の P/ρ/(γ-1) は TP では e(T) と不一致になるため使わない。
            const flow_float e0c = roe_c[ic0]/ro[ic0]
                                 - 0.5*(Ux[ic0]*Ux[ic0]+Uy[ic0]*Uy[ic0]+Uz[ic0]*Uz[ic0]);
            const flow_float e1c = roe_c[ic1]/ro[ic1]
                                 - 0.5*(Ux[ic1]*Ux[ic1]+Uy[ic1]*Uy[ic1]+Uz[ic1]*Uz[ic1]);
            Itilde = Ctilde*0.5*(e0c + e1c);
        } else {
            Itilde = Ctilde*0.5*(Ps[ic0]/ro[ic0] +Ps[ic1]/ro[ic1])/(ga-1.0);   // CPG (ビット不変)
        }
        // free-stream 保存: 運動量圧力項は基準静圧 d_pRef を面ごとに差し引いて組む
        // (定数ゲージシフト: Σs=0 解析的に解不変。float32 の p·s 桁落ちによる偽運動量源を除去。
        //  SLAU と同じ処方 = plans/active/convection-freestream-preserving-flux.md。
        //  エネルギー圧力仕事 Ptilde は -pRef∇·u で方程式が変わるため差分しない。pRef=0 でビット不変)
        flow_float Gtildex = 0.5*((Ps[ic0]-d_pRef)+(Ps[ic1]-d_pRef))*nx;
        flow_float Gtildey = 0.5*((Ps[ic0]-d_pRef)+(Ps[ic1]-d_pRef))*ny;
        flow_float Gtildez = 0.5*((Ps[ic0]-d_pRef)+(Ps[ic1]-d_pRef))*nz;
        flow_float Ptilde  = 0.5*( (Ux[ic0]*Ps[ic1] + Ux[ic1]*Ps[ic0])*nx
                                  +(Uy[ic0]*Ps[ic1] + Uy[ic1]*Ps[ic0])*ny
                                  +(Uz[ic0]*Ps[ic1] + Uz[ic1]*Ps[ic0])*nz );

        res_ro_temp   = Ctilde*sss;
        res_roe_temp  = (Ktilde + Itilde + Ptilde)*sss;
        res_roUx_temp = (Mtildex + Gtildex)*sss;
        res_roUy_temp = (Mtildey + Gtildey)*sss;
        res_roUz_temp = (Mtildez + Gtildez)*sss;
        } else {
        // ---- advGauge: 差分形因数分解 (小さい量を先に作る; 一様場で各因子がビット単位ゼロ) ----
        const flow_float d0x = Ux[ic0]-d_uRefX, d0y = Uy[ic0]-d_uRefY, d0z = Uz[ic0]-d_uRefZ;
        const flow_float d1x = Ux[ic1]-d_uRefX, d1y = Uy[ic1]-d_uRefY, d1z = Uz[ic1]-d_uRefZ;
        const flow_float dro = 0.5*((ro[ic0]-d_roRef)+(ro[ic1]-d_roRef));
        const flow_float dux = 0.5*(d0x+d1x), duy = 0.5*(d0y+d1y), duz = 0.5*(d0z+d1z);
        const flow_float uxb = 0.5*(Ux[ic0]+Ux[ic1]);
        const flow_float uyb = 0.5*(Uy[ic0]+Uy[ic1]);
        const flow_float uzb = 0.5*(Uz[ic0]+Uz[ic1]);
        const flow_float Unb = uxb*nx + uyb*ny + uzb*nz;
        const flow_float dUn = dux*nx + duy*ny + duz*nz;
        const flow_float UnR = d_uRefX*nx + d_uRefY*ny + d_uRefZ*nz;
        const flow_float Cinf = d_roRef*UnR;
        const flow_float dC   = dro*Unb + d_roRef*dUn;              // = Ctilde − C∞
        // 圧力項 (p−pRef) は従来 Gtilde と同一
        const flow_float dpr0 = Ps[ic0]-d_pRef, dpr1 = Ps[ic1]-d_pRef;
        const flow_float Gp = 0.5*(dpr0+dpr1);
        // エネルギー: K̃, Ĩ, P̃ を各々 (差分)×(平均)+(参照)×(差分) で
        const flow_float kb = 0.5*(Ux[ic0]*Ux[ic1] +Uy[ic0]*Uy[ic1] +Uz[ic0]*Uz[ic1]);
        const flow_float kdiff = 0.5*( d0x*Ux[ic1] + d0y*Uy[ic1] + d0z*Uz[ic1]
                                      + d_uRefX*d1x + d_uRefY*d1y + d_uRefZ*d1z );   // = k̄ − ½|u∞|²
        const flow_float eRefG = d_pRef/d_roRef;                    // (γ−1)e∞ (同一 float 除算でビット整合)
        const flow_float pr0 = Ps[ic0]/ro[ic0], pr1 = Ps[ic1]/ro[ic1];
        const flow_float ebar = 0.5*(pr0+pr1)/(ga-1.0);
        const flow_float de   = 0.5*((pr0-eRefG)+(pr1-eRefG))/(ga-1.0);
        const flow_float dP = 0.5*( (d0x*Ps[ic1] + d_uRefX*dpr1 + d1x*Ps[ic0] + d_uRefX*dpr0)*nx
                                   +(d0y*Ps[ic1] + d_uRefY*dpr1 + d1y*Ps[ic0] + d_uRefY*dpr0)*ny
                                   +(d0z*Ps[ic1] + d_uRefZ*dpr1 + d1z*Ps[ic0] + d_uRefZ*dpr0)*nz );

        CinfS = Cinf*sss;
        res_ro_temp   = dC*sss;
        res_roUx_temp = (dC*uxb + Cinf*dux + Gp*nx)*sss;
        res_roUy_temp = (dC*uyb + Cinf*duy + Gp*ny)*sss;
        res_roUz_temp = (dC*uzb + Cinf*duz + Gp*nz)*sss;
        res_roe_temp  = (dC*kb + Cinf*kdiff + dC*ebar + Cinf*de + dP)*sss;
        }

        // ---- opt-in entropy-stable 散逸レイヤ (keepDissType 1: scalar / 2: matrix) ----
        //   1 (scalar): F -= 0.5*σ*λ'*ΔU。LLF 型は Δw·ΔU>=0 で ES。全成分同一 λ' なので
        //     市松と一緒に渦も食う (TGV KE cost 2.7%@σ=0.05)。
        //   2 (matrix): F -= 0.5*σ*R|Λ'|S Rᵀ Δw (entropy-scaled Roe 型, Chandrashekar/Barth)。
        //     w=∂η/∂U (η=-ρs/(γ-1)), H=∂U/∂w=RSRᵀ (S: 音響 ρ/(2γ)・エントロピー ρ(γ-1)/γ・
        //     せん断 p — tools/verify_entropy_scaling.py で数値検証済)。
        //     |Λ'| は音響のみ |Un|+c'、せん断/エントロピーは |Un| → 渦を食わず市松 (音響) を減衰。
        //   λ'/c': lowMachPrecond>=1 なら c'=lowMachCprime (低マッハで c'→O(|u|)、M>=1 で c に復帰)。
        //   過大 σ は解像乱流を殺すため TGV (L2) で較正する。
        //   f_d 駆動 σ ブレンド (fd_shield!=nullptr): σ_f=max(σ_min,(1-f̃_d)·σ_max)。RANS 帯 (f̃_d≈0)
        //   でフル ES 散逸 (他コードの風上相当)、LES 域はフロアのみ (Travin/SU2 FD/UZEN 型ゾーニング)。
        //   halo/ghost セル (ic1>=nCells) は fd_shield 未計算のため owner 値で代表する。
        flow_float sigmaF = keepDissCoeff;
        if (fd_shield != nullptr) {
            const flow_float fd0 = fd_shield[ic0];
            const flow_float fd1 = (ic1 < nCells) ? fd_shield[ic1] : fd0;
            const flow_float fdF = 0.5*(fd0 + fd1);
            // ★ fd_shield の向きは DESmode で逆 (variables.hpp): DDES=f_d (1=LES) / IDDES=f̃_d (1=RANS)
            const flow_float ransFrac = (fdLesOne == 1) ? (1.0 - fdF) : fdF;
            sigmaF = fmax(keepDissCoeff, ransFrac*keepDissCoeffMax);
        }
        if (keepDissType == 1) {
            // 保存量ジャンプ (primitives から構成; CPG: roe = P/(γ-1) + ½ρ|u|²)
            const flow_float dro   = ro[ic1] - ro[ic0];
            const flow_float droUx = ro[ic1]*Ux[ic1] - ro[ic0]*Ux[ic0];
            const flow_float droUy = ro[ic1]*Uy[ic1] - ro[ic0]*Uy[ic0];
            const flow_float droUz = ro[ic1]*Uz[ic1] - ro[ic0]*Uz[ic0];
            // TP では保存量 roe のジャンプそのもの (EOS 非依存)、CPG は従来式 (ビット不変)
            flow_float droe;
            if (thermalMethod == 2) {
                droe = roe_c[ic1] - roe_c[ic0];
            } else {
                const flow_float roe0  = Ps[ic0]/(ga-1.0)
                                       + 0.5*ro[ic0]*(Ux[ic0]*Ux[ic0]+Uy[ic0]*Uy[ic0]+Uz[ic0]*Uz[ic0]);
                const flow_float roe1  = Ps[ic1]/(ga-1.0)
                                       + 0.5*ro[ic1]*(Ux[ic1]*Ux[ic1]+Uy[ic1]*Uy[ic1]+Uz[ic1]*Uz[ic1]);
                droe = roe1 - roe0;
            }

            // 面代表の波速 λ' = |Un| + c(') (面平均状態で評価)
            const flow_float roF = 0.5*(ro[ic0]+ro[ic1]);
            const flow_float PsF = 0.5*(Ps[ic0]+Ps[ic1]);
            const flow_float uxF = 0.5*(Ux[ic0]+Ux[ic1]);
            const flow_float uyF = 0.5*(Uy[ic0]+Uy[ic1]);
            const flow_float uzF = 0.5*(Uz[ic0]+Uz[ic1]);
            const flow_float Un  = uxF*nx + uyF*ny + uzF*nz;
            const flow_float c   = (thermalMethod == 2)
                                 ? 0.5*(st.sonic[ic0]+st.sonic[ic1])   // TP: dependent vars の凍結音速
                                 : sqrt(ga*PsF/roF);                   // CPG (ビット不変)
            const flow_float cd  = (keepDissCprime == 1)
                                 ? lowMachCprime(c, sqrt(uxF*uxF+uyF*uyF+uzF*uzF), Un, precondEps)
                                 : c;
            const flow_float coef = 0.5*sigmaF*(fabs(Un) + cd)*sss;

            res_ro_temp   -= coef*dro;
            res_roUx_temp -= coef*droUx;
            res_roUy_temp -= coef*droUy;
            res_roUz_temp -= coef*droUz;
            res_roe_temp  -= coef*droe;
        }
        else if (keepDissType == 2 && thermalMethod != 2) {
            // ---- matrix ES 散逸 (CPG): D = R|Λ'|S Rᵀ Δw ----
            // keepDissJump==1 (再構成後ジャンプ): Δw を面へ線形再構成した L/R 状態
            //   q_L = q_i + g_i·(x_f−x_i), q_R = q_j + g_j·(x_f−x_j) (κ=0・リミタ無し) から組む。
            //   純 2Δ モード (市松) は中心勾配が厳密ゼロ → Δ 不変でフル減衰のまま、滑らかな場は
            //   Δ が O(h)→O(h³) に落ち解像スケールのドレインを桁減 (実測: 市松 1.00 / TGV 層流期
            //   0.006 / 遷移期 0.16 / ピーク 0.45。plans/active/convection-keep-diss-recon-jump.md)。
            //   ghost 面 (cell の ip>=nNormalPlanes) は ghost 勾配/重心が無効なので生ジャンプへ
            //   フォールバック。node の主ループは内部双対面のみで勾配 gather 済 → 常に再構成可。
            //   再構成で ρ/p が非正になる面も生ジャンプへ (リミタ無しの正値性ガード)。
            //   面平均量 (roF, c, Ht 等) は従来どおりセル生値 (ビット不変性は keepDissJump=0 で担保)。
            //   keepDissJump==2 (sign-property, TeCNO 型): 再構成ジャンプの特性射影 rd_k を
            //   生ジャンプ射影と minmod し、各波で rd_used·rd_raw >= 0 を構造的に保証
            //   → エントロピー散逸性 Σ S_k λ_k rd_used_k rd_raw_k >= 0 が証明付きで復活
            //   (Fjordholm-Mishra-Tadmor の sign property の射影クリップ版)。
            //   事前測定: 符号反転は (面,波) の 27-37% だが負寄与は総散逸の 1-2% で、
            //   minmod 後の散逸総量は jump=1 とほぼ同一 (市松 0.996 / 遷移期 0.32)。
            flow_float ro0j = ro[ic0], ro1j = ro[ic1];
            flow_float ux0j = Ux[ic0], uy0j = Uy[ic0], uz0j = Uz[ic0];
            flow_float ux1j = Ux[ic1], uy1j = Uy[ic1], uz1j = Uz[ic1];
            flow_float p0j  = Ps[ic0], p1j  = Ps[ic1];
            bool cbRecon = false;   // cb 補正の有効面フラグ (内部面かつ両側再構成成功のみ)
            if (keepDissJump >= 1 && ip < geom.nNormalPlanes) {
                // pc − cc[ic0] / pc − cc[ic1] は double の値の位置で引いて 1 回だけ丸めた r0 / r1 を読む (段 ④、plan
                // architecture-float-state-double-geometry §4.2a)
                const geom_float d0x = geom.gr0_x[ip], d0y = geom.gr0_y[ip], d0z = geom.gr0_z[ip];
                const geom_float d1x = geom.gr1_x[ip], d1y = geom.gr1_y[ip], d1z = geom.gr1_z[ip];
                const flow_float roL = ro[ic0] + grd.drodx[ic0]*d0x + grd.drody[ic0]*d0y + grd.drodz[ic0]*d0z;
                const flow_float roR = ro[ic1] + grd.drodx[ic1]*d1x + grd.drody[ic1]*d1y + grd.drodz[ic1]*d1z;
                const flow_float pL  = Ps[ic0] + grd.dPdx[ic0]*d0x + grd.dPdy[ic0]*d0y + grd.dPdz[ic0]*d0z;
                const flow_float pR  = Ps[ic1] + grd.dPdx[ic1]*d1x + grd.dPdy[ic1]*d1y + grd.dPdz[ic1]*d1z;
                if (roL > 0.0 && roR > 0.0 && pL > 0.0 && pR > 0.0) {
                    cbRecon = true;
                    ro0j = roL; ro1j = roR; p0j = pL; p1j = pR;
                    ux0j = Ux[ic0] + grd.dUxdx[ic0]*d0x + grd.dUxdy[ic0]*d0y + grd.dUxdz[ic0]*d0z;
                    uy0j = Uy[ic0] + grd.dUydx[ic0]*d0x + grd.dUydy[ic0]*d0y + grd.dUydz[ic0]*d0z;
                    uz0j = Uz[ic0] + grd.dUzdx[ic0]*d0x + grd.dUzdy[ic0]*d0y + grd.dUzdz[ic0]*d0z;
                    ux1j = Ux[ic1] + grd.dUxdx[ic1]*d1x + grd.dUxdy[ic1]*d1y + grd.dUxdz[ic1]*d1z;
                    uy1j = Uy[ic1] + grd.dUydx[ic1]*d1x + grd.dUydy[ic1]*d1y + grd.dUydz[ic1]*d1z;
                    uz1j = Uz[ic1] + grd.dUzdx[ic1]*d1x + grd.dUzdy[ic1]*d1y + grd.dUzdz[ic1]*d1z;
                }
            }
            // Δw (エントロピー変数ジャンプ): w=[(γ-s)/(γ-1)-β|u|², 2βu, -2β], β=ρ/(2p), s=ln p-γ ln ρ
            const flow_float b0 = ro0j/(2.0*p0j);
            const flow_float b1 = ro1j/(2.0*p1j);
            const flow_float q0 = ux0j*ux0j+uy0j*uy0j+uz0j*uz0j;
            const flow_float q1 = ux1j*ux1j+uy1j*uy1j+uz1j*uz1j;
            const flow_float ds = log(p1j/p0j) - ga*log(ro1j/ro0j);
            const flow_float dw0 = -ds/(ga-1.0) - (b1*q1 - b0*q0);
            const flow_float dw1 = 2.0*(b1*ux1j - b0*ux0j);
            const flow_float dw2 = 2.0*(b1*uy1j - b0*uy0j);
            const flow_float dw3 = 2.0*(b1*uz1j - b0*uz0j);
            const flow_float dw4 = -2.0*(b1 - b0);

            // 面平均状態と固有系 (python 検証済の正規化: r=[1,u∓cn,Ht∓c·Un] 等)
            const flow_float roF = 0.5*(ro[ic0]+ro[ic1]);
            const flow_float PsF = 0.5*(Ps[ic0]+Ps[ic1]);
            const flow_float uxF = 0.5*(Ux[ic0]+Ux[ic1]);
            const flow_float uyF = 0.5*(Uy[ic0]+Uy[ic1]);
            const flow_float uzF = 0.5*(Uz[ic0]+Uz[ic1]);
            const flow_float qF  = uxF*uxF+uyF*uyF+uzF*uzF;
            const flow_float Un  = uxF*nx + uyF*ny + uzF*nz;
            const flow_float c   = sqrt(ga*PsF/roF);
            const flow_float Ht  = ga*PsF/((ga-1.0)*roF) + 0.5*qF;

            // 接線基底 (n に対しロバストに構築)
            flow_float ax, ay, az;
            if (fabs(nx) < 0.9) { ax=1.0; ay=0.0; az=0.0; } else { ax=0.0; ay=1.0; az=0.0; }
            flow_float t1x = ny*az - nz*ay, t1y = nz*ax - nx*az, t1z = nx*ay - ny*ax;
            const flow_float t1n = sqrt(t1x*t1x+t1y*t1y+t1z*t1z);
            t1x/=t1n; t1y/=t1n; t1z/=t1n;
            const flow_float t2x = ny*t1z - nz*t1y, t2y = nz*t1x - nx*t1z, t2z = nx*t1y - ny*t1x;

            // 固有値: 音響のみ c(') 込み、せん断/エントロピーは |Un| (渦を守る)
            const flow_float cd  = (keepDissCprime == 1)
                                 ? lowMachCprime(c, sqrt(qF), Un, precondEps)
                                 : c;
            const flow_float lamA = fabs(Un) + cd;   // 音響 (±両波共通の安全上界)
            const flow_float lamU = fabs(Un);        // エントロピー・せん断

            // スケーリング S (Barth): 音響 ρ/(2γ), エントロピー ρ(γ-1)/γ, せん断 p
            const flow_float sA = roF/(2.0*ga);
            const flow_float sE = roF*(ga-1.0)/ga;
            const flow_float sS = PsF;

            // z_k = S_k |λ_k| (r_k · Δw)
            const flow_float ut1 = uxF*t1x + uyF*t1y + uzF*t1z;
            const flow_float ut2 = uxF*t2x + uyF*t2y + uzF*t2z;
            flow_float rd1 = dw0 + (uxF-c*nx)*dw1 + (uyF-c*ny)*dw2 + (uzF-c*nz)*dw3 + (Ht-c*Un)*dw4; // r(un-c)·Δw
            flow_float rd2 = dw0 + uxF*dw1 + uyF*dw2 + uzF*dw3 + 0.5*qF*dw4;                          // r(entropy)·Δw
            flow_float rd3 = t1x*dw1 + t1y*dw2 + t1z*dw3 + ut1*dw4;                                   // r(shear1)·Δw
            flow_float rd4 = t2x*dw1 + t2y*dw2 + t2z*dw3 + ut2*dw4;                                   // r(shear2)·Δw
            flow_float rd5 = dw0 + (uxF+c*nx)*dw1 + (uyF+c*ny)*dw2 + (uzF+c*nz)*dw3 + (Ht+c*Un)*dw4; // r(un+c)·Δw
            // opBlend 標準 Roe 側用の「生ジャンプ」音響射影。既定は rd と同一 (=従来挙動)。
            // jump==2 のとき下で rr1/rr5 (真の生ジャンプ) に置換される。
            flow_float rdAraw1 = rd1, rdAraw5 = rd5;
            if (keepDissJump == 2) {
                // sign-property クリップ: 生ジャンプ Δw_raw の特性射影と成分ごと minmod。
                // 符号一致なら小さい方 (通常は再構成側=高次の利得維持)、反転なら 0 (ES 保証)。
                // 再構成が発動しなかった面 (ghost/正値性 fallback) は dw==dw_raw で恒等。
                const flow_float br0 = ro[ic0]/(2.0*Ps[ic0]);
                const flow_float br1 = ro[ic1]/(2.0*Ps[ic1]);
                const flow_float qr0 = Ux[ic0]*Ux[ic0]+Uy[ic0]*Uy[ic0]+Uz[ic0]*Uz[ic0];
                const flow_float qr1 = Ux[ic1]*Ux[ic1]+Uy[ic1]*Uy[ic1]+Uz[ic1]*Uz[ic1];
                const flow_float dsr = log(Ps[ic1]/Ps[ic0]) - ga*log(ro[ic1]/ro[ic0]);
                const flow_float dwr0 = -dsr/(ga-1.0) - (br1*qr1 - br0*qr0);
                const flow_float dwr1 = 2.0*(br1*Ux[ic1] - br0*Ux[ic0]);
                const flow_float dwr2 = 2.0*(br1*Uy[ic1] - br0*Uy[ic0]);
                const flow_float dwr3 = 2.0*(br1*Uz[ic1] - br0*Uz[ic0]);
                const flow_float dwr4 = -2.0*(br1 - br0);
                const flow_float rr1 = dwr0 + (uxF-c*nx)*dwr1 + (uyF-c*ny)*dwr2 + (uzF-c*nz)*dwr3 + (Ht-c*Un)*dwr4;
                const flow_float rr2 = dwr0 + uxF*dwr1 + uyF*dwr2 + uzF*dwr3 + 0.5*qF*dwr4;
                const flow_float rr3 = t1x*dwr1 + t1y*dwr2 + t1z*dwr3 + ut1*dwr4;
                const flow_float rr4 = t2x*dwr1 + t2y*dwr2 + t2z*dwr3 + ut2*dwr4;
                const flow_float rr5 = dwr0 + (uxF+c*nx)*dwr1 + (uyF+c*ny)*dwr2 + (uzF+c*nz)*dwr3 + (Ht+c*Un)*dwr4;
                rd1 = (rd1*rr1 <= 0.0) ? 0.0 : (fabs(rd1) < fabs(rr1) ? rd1 : rr1);
                rd2 = (rd2*rr2 <= 0.0) ? 0.0 : (fabs(rd2) < fabs(rr2) ? rd2 : rr2);
                rd3 = (rd3*rr3 <= 0.0) ? 0.0 : (fabs(rd3) < fabs(rr3) ? rd3 : rr3);
                rd4 = (rd4*rr4 <= 0.0) ? 0.0 : (fabs(rd4) < fabs(rr4) ? rd4 : rr4);
                rd5 = (rd5*rr5 <= 0.0) ? 0.0 : (fabs(rd5) < fabs(rr5) ? rd5 : rr5);
                rdAraw1 = rr1; rdAraw5 = rr5;   // opBlend 標準 Roe 側の生ジャンプ (下記参照)
            }
            flow_float z1, z5;
            // fdblend×precond の演算子ブレンド (X4): 実効散逸 = σ_min·D_precond + (σ_f−σ_min)·D_std。
            //   LES 域 (σ_f=σ_min) は純 precond (較正済み市松キラー)、壁帯 (σ_f→1) は precond 寄与を
            //   σ_min に固定したまま標準 Roe を増やす → Δp 増強 (∝c²/Ur)×σ→1 の爆発を回避しつつ
            //   Δp 減衰は常時確保。coef=0.5·σ_f·sss が外掛けなので比率 r=σ_min/σ_f で内分する。
            const bool opBlend = (keepDissPrecond == 1) && (fd_shield != nullptr);
            const flow_float rOp = opBlend ? (keepDissCoeff / sigmaF) : 1.0;
            if (keepDissPrecond == 1) {
                // ---- Turkel 前処理音響散逸 (plans/active/convection-keep-diss-lowmach-precond.md §3) ----
                // α∓ = S_A rd∓ (clip/再構成適用済) を特性成分 (Δp, ΔUn) に写像し、前処理 2×2
                //   D_p = Γ|Γ⁻¹A2| (閉形 |M2| = φ1 M2 + φ2 I, 固有ベクトル不要) を掛けて z∓ に戻す。
                // 漸近: Δp 散逸 ∝ c²/Ur (連続式への圧力 Laplacian = Rhie-Chow 相当の市松キラー)、
                //   ΔUn 散逸 ∝ Ur (Guillard-Viozat の低マッハ過散逸を解消)。β=1 (M>=1) で標準 Roe |A2| に復帰。
                // sym(K) は全 (M,Un) で正定値 → ES 性維持 (tools/verify_precond_dissipation.py 全 6 検証 PASS,
                //   完全 Weiss-Smith 5×5 との一致 6e-4)。keepDissCprime は本枝では不使用 (前処理が上位互換)。
                const flow_float velMag = sqrt(qF);
                const flow_float Ur   = lowMachUr(c, velMag, precondEps);
                const flow_float beta = (Ur/c)*(Ur/c);
                const flow_float up   = 0.5*(1.0+beta)*Un;
                const flow_float cp   = lowMachCprime(c, velMag, Un, precondEps);
                const flow_float lp = up + cp, lm = up - cp;
                const flow_float phi1 = (fabs(lp) - fabs(lm))/(2.0*cp);
                const flow_float phi2 = (lp*fabs(lm) - lm*fabs(lp))/(2.0*cp);
                const flow_float d11 = phi1*Un + phi2/beta;
                const flow_float d12 = phi1*roF*c*c;
                const flow_float d21 = phi1/roF;
                const flow_float d22 = phi1*Un + phi2;
                const flow_float am = sA*rd1, ap = sA*rd5;
                const flow_float dpc = c*c*(ap + am);        // Δp (特性成分)
                const flow_float dun = c*(ap - am)/roF;      // ΔUn (特性成分)
                const flow_float fp = d11*dpc + d12*dun;
                const flow_float fu = d21*dpc + d22*dun;
                z1 = 0.5*(fp - roF*c*fu)/(c*c);
                z5 = 0.5*(fp + roF*c*fu)/(c*c);
                if (opBlend) {
                    // 内分: r·precond + (1−r)·標準 Roe (r=σ_min/σ_f, LES 域で r=1=純 precond)。
                    // 標準 Roe 側は **フル c** (c' でなく) を使う: 増分が入るのは RANS/grey 帯のみで、
                    // そこは渦保護不要、かつ大波長音響モードの減衰には音響スケール ~c が必須
                    // (P3/P4 減衰試験: c' のままだと減衰半分。run_diag_pdamp_* 2026-07-22)。
                    // さらに keepDissOpBlendRaw=1 なら標準 Roe 側は**生ジャンプ** (rdAraw, jump2 の
                    // 再構成フィルタ非適用) を使う: RANS/grey 帯で 2-4Δ を「解像スケール」として
                    // 保護する理由はなく (DES 定石=RANS 帯は完全風上相当)、jump2 フィルタが
                    // 掛かったままだと σ→1 にしても 2-3Δ 定在モードへの実効ジャンプが中央値 6% に
                    // 落ちて散逸ゾーニングが機能しない (case/39 丘頂鋸歯, run_0015-0018 掃引)。
                    // 生ジャンプは一様場で厳密 0 → free-stream 安全。
                    const flow_float lamAfull = fabs(Un) + c;
                    const flow_float rdA1 = (keepDissOpBlendRaw == 1) ? rdAraw1 : rd1;
                    const flow_float rdA5 = (keepDissOpBlendRaw == 1) ? rdAraw5 : rd5;
                    z1 = rOp*z1 + (1.0-rOp)*sA*lamAfull*rdA1;
                    z5 = rOp*z5 + (1.0-rOp)*sA*lamAfull*rdA5;
                }
            } else {
                z1 = sA*lamA*rd1;
                z5 = sA*lamA*rd5;
            }
            const flow_float z2 = sE*lamU*rd2;
            const flow_float z3 = sS*lamU*rd3;
            const flow_float z4 = sS*lamU*rd4;

            // D = Σ z_k r_k、flux -= 0.5σ D S
            const flow_float coef = 0.5*sigmaF*sss;
            const flow_float Dro   = z1 + z2 + z5;
            const flow_float DroUx = z1*(uxF-c*nx) + z2*uxF + z3*t1x + z4*t2x + z5*(uxF+c*nx);
            const flow_float DroUy = z1*(uyF-c*ny) + z2*uyF + z3*t1y + z4*t2y + z5*(uyF+c*ny);
            const flow_float DroUz = z1*(uzF-c*nz) + z2*uzF + z3*t1z + z4*t2z + z5*(uzF+c*nz);
            const flow_float Droe  = z1*(Ht-c*Un) + z2*0.5*qF + z3*ut1 + z4*ut2 + z5*(Ht+c*Un);

            res_ro_temp   -= coef*Dro;
            res_roUx_temp -= coef*DroUx;
            res_roUy_temp -= coef*DroUy;
            res_roUz_temp -= coef*DroUz;
            res_roe_temp  -= coef*Droe;

            // ---- 高周波圧力欠陥駆動 mass-flux 補正 (Rhie–Chow 型市松キラー, opt-in) ----
            // plans/active/convection-keep-cb-pressure-correction.md。
            // δp^HF = pR−pL (再構成不整合): 線形圧力場で厳密ゼロ・2Δ 市松で生 Δp → 市松狙い撃ち。
            // 有効面 = 内部面かつ両側再構成成功のみ (生 Δp フォールバック禁止: 滑らかな物理勾配を
            // 拡散させない)。σ_min×precond の Δp 項への加算で、σ/f̃_d/opBlend から独立。
            // sign gate: rank-1 flux −k δp^HF g (g=[1,ū,H̄t]) の entropy 寄与は −k δp^HF (Δw·g)
            // なので δp^HF·r_f>0 の面のみ通し、面ごとに散逸性を構造保証
            // (tools/verify_cb_pressure_correction.py: 市松素通し 100%・乱流 10% ジャンプで gate 落ち 0.2%)。
            // 全保存量へ ū,H̄t で整合配布 → massflux (散逸込み) 経由でスカラー/種/乱流も自動整合。
            if (keepDissCbCoeff > 0.0 && cbRecon) {
                // δp^HF は**差分形**で組む: (Ps1−Ps0) + (g1·d1 − g0·d0)。絶対圧力の再構成
                // pR−pL (数学的に同値) は float32 で p~1e5 の量子化ノイズを拾う。
                // pc − cc は r0 / r1 を読む (段 ④、上と同じ)
                const geom_float e0x = geom.gr0_x[ip], e0y = geom.gr0_y[ip], e0z = geom.gr0_z[ip];
                const geom_float e1x = geom.gr1_x[ip], e1y = geom.gr1_y[ip], e1z = geom.gr1_z[ip];
                const flow_float dpRaw = Ps[ic1] - Ps[ic0];
                const flow_float dpRec = dpRaw
                                       + (grd.dPdx[ic1]*e1x + grd.dPdy[ic1]*e1y + grd.dPdz[ic1]*e1z)
                                       - (grd.dPdx[ic0]*e0x + grd.dPdy[ic0]*e0y + grd.dPdz[ic0]*e0z);
                // sign-property minmod (jump2 の rd クリップと同原理): 生 Δp と符号一致なら小さい方、
                // 反転なら 0。一様場では dpRaw が**厳密 0** → GG 勾配のメトリック閉合ノイズ
                // (g·d ~ 0.1 Pa 級、1/Ur 増幅で free-stream を 5 step で破壊: case/33 run_0021/0022)
                // を構造的に完全消去。2Δ 市松は中心勾配 0 → dpRec=dpRaw でフル振幅のまま。
                const flow_float dpHF = (dpRec*dpRaw <= 0.0) ? 0.0
                                      : (fabs(dpRec) < fabs(dpRaw) ? dpRec : dpRaw);
                const flow_float rf = dw0 + uxF*dw1 + uyF*dw2 + uzF*dw3 + Ht*dw4;  // Δw·g
                if (dpHF*rf > 0.0) {
                    const flow_float UrCb = lowMachUr(c, sqrt(qF), keepDissCbEps);
                    const flow_float mcb  = 0.5*keepDissCbCoeff*sss*dpHF/UrCb;     // δṁ (符号は高圧→低圧)
                    res_ro_temp   -= mcb;
                    res_roUx_temp -= mcb*uxF;
                    res_roUy_temp -= mcb*uyF;
                    res_roUz_temp -= mcb*uzF;
                    res_roe_temp  -= mcb*Ht;
                }
            }
        }
        else if (keepDissType == 2) {
            // ---- matrix ES 散逸 (TP, Step 3/4): D = R|Λ'|S Rᵀ Δw ----
            // 次元付きエントロピー η=-ρs_mix, w=[(g-½|u|²)/T, u/T, -1/T] (Chalot-Hughes-Shakib),
            // S: 音響 ρ/(2γR_mix)・エントロピー ρ/cp_mix(T,Y)・せん断 ρT、エントロピー波 r_E=½q+e-cv·T
            // (tools/verify_entropy_scaling_tp.py + 凍結組成混合の数値検証で H=RSRᵀ・Ar=λr・SPD 確認済)。
            // 多成分 (Step 4): 混合物性 (R_mix/cp_mix/s0_mix) + 混合エントロピー -ΣY_kR_k ln X_k
            //   (Y ln X → 0 でゼロ濃度種でも正則)。組成一様なら凍結組成の TP 1 気体に厳密縮退。
            //   組成ジャンプ面では全系の厳密 ES は主張しない (Q は SPD なので無条件に散逸的=安定)。
            //   種輸送は散逸込み massflux 経由で連続式と整合 (種方程式自体への行列散逸は将来課題)。
            // thermo 評価は double (thermoHrefTemp datum 焼き込み係数をそのまま使用 → w も自動整合)。
            const int ns = (nSpecies > 1) ? nSpecies : 1;
            double Y0[THERMO_MAX_SPECIES], Y1[THERMO_MAX_SPECIES], YF[THERMO_MAX_SPECIES];
            if (nSpecies > 1) {
                // ★ Y はカーネル内で正規化 (Σ_k Y_k = 1 を強制)。ρY と ρ は別カーネル/atomicAdd で
                //   更新されるため Σ ρY_k ≠ ρ の共通モードノイズ (~1e-7) が乗り、これが s⁰ (×~10) と
                //   ln X (対数微分特異) で増幅されて w ノイズ→弱減衰モードへ注入され、市松減衰が
                //   プラトー化する (Y=[1,0] 縮退 bisect で確定)。正規化で共通モードを除去する。
                double sum0 = 0.0, sum1 = 0.0;
                for (int k = 0; k < ns; k++) {
                    Y0[k] = fmax((double)roY[k][ic0], 0.0);
                    Y1[k] = fmax((double)roY[k][ic1], 0.0);
                    sum0 += Y0[k]; sum1 += Y1[k];
                }
                const double inv0 = (sum0 > 1e-300) ? 1.0/sum0 : 0.0;
                const double inv1 = (sum1 > 1e-300) ? 1.0/sum1 : 0.0;
                for (int k = 0; k < ns; k++) {
                    Y0[k] *= inv0; Y1[k] *= inv1;
                    YF[k] = 0.5*(Y0[k]+Y1[k]);
                }
                if (sum0 <= 1e-300) { Y0[0] = 1.0; YF[0] = 0.5*(Y0[0]+Y1[0]); }
                if (sum1 <= 1e-300) { Y1[0] = 1.0; YF[0] = 0.5*(Y0[0]+Y1[0]); }
            } else { Y0[0] = 1.0; Y1[0] = 1.0; YF[0] = 1.0; }
            const double R0g = thermo_R_mix(sp, ns, Y0);
            const double R1g = thermo_R_mix(sp, ns, Y1);
            const double q0d = (double)Ux[ic0]*Ux[ic0]+(double)Uy[ic0]*Uy[ic0]+(double)Uz[ic0]*Uz[ic0];
            const double q1d = (double)Ux[ic1]*Ux[ic1]+(double)Uy[ic1]*Uy[ic1]+(double)Uz[ic1]*Uz[ic1];
            const double T0d = (double)Ps[ic0]/((double)ro[ic0]*R0g);
            const double T1d = (double)Ps[ic1]/((double)ro[ic1]*R1g);
            const double PsFd = 0.5*((double)Ps[ic0]+(double)Ps[ic1]);
            // 混合エントロピー項 -Σ Y_k R_k ln X_k (単成分では 0)
            double mixent0 = 0.0, mixent1 = 0.0;
            if (nSpecies > 1) {
                double mol0 = 0.0, mol1 = 0.0;
                for (int k = 0; k < ns; k++) { mol0 += Y0[k]/sp[k].MW; mol1 += Y1[k]/sp[k].MW; }
                for (int k = 0; k < ns; k++) {
                    const double Rk = THERMO_RU/sp[k].MW;
                    if (Y0[k] > 0.0) mixent0 -= Y0[k]*Rk*log(fmax(Y0[k]/sp[k].MW/mol0, 1e-300));
                    if (Y1[k] > 0.0) mixent1 -= Y1[k]*Rk*log(fmax(Y1[k]/sp[k].MW/mol1, 1e-300));
                }
            }
            // s_i: log 引数 ~1 で桁落ち回避のため p_i/PsF で評価し、多成分では NASA-9 標準状態
            //   p°=THERMO_P_STD への datum 補正 −R_i·ln(PsF/p°) を加える。面ローカル pref=PsF のままだと
            //   組成ジャンプ面 (R0≠R1) で (R1−R0)·ln(PsF) が Δw に残り、散逸量が面参照圧という任意の
            //   選択に依存してしまう (単成分は R0=R1 で解析相殺するため補正 0.0 加算=ビット不変)。
            const double lnPstd = (nSpecies > 1) ? log(PsFd/THERMO_P_STD) : 0.0;
            const double s0d = thermo_s0_mix(sp, ns, Y0, T0d) + mixent0 - R0g*(log((double)Ps[ic0]/PsFd) + lnPstd);
            const double s1d = thermo_s0_mix(sp, ns, Y1, T1d) + mixent1 - R1g*(log((double)Ps[ic1]/PsFd) + lnPstd);
            const double g0d = thermo_h_mix(sp, ns, Y0, T0d) - T0d*s0d;
            const double g1d = thermo_h_mix(sp, ns, Y1, T1d) - T1d*s1d;
            const flow_float dw0 = (flow_float)((g1d-0.5*q1d)/T1d - (g0d-0.5*q0d)/T0d);
            const flow_float dw1 = (flow_float)((double)Ux[ic1]/T1d - (double)Ux[ic0]/T0d);
            const flow_float dw2 = (flow_float)((double)Uy[ic1]/T1d - (double)Uy[ic0]/T0d);
            const flow_float dw3 = (flow_float)((double)Uz[ic1]/T1d - (double)Uz[ic0]/T0d);
            const flow_float dw4 = (flow_float)(-(1.0/T1d - 1.0/T0d));

            // 面平均状態 (γ,cp,cv は面平均温度で凍結評価)
            const flow_float roF = 0.5*(ro[ic0]+ro[ic1]);
            const flow_float uxF = 0.5*(Ux[ic0]+Ux[ic1]);
            const flow_float uyF = 0.5*(Uy[ic0]+Uy[ic1]);
            const flow_float uzF = 0.5*(Uz[ic0]+Uz[ic1]);
            const flow_float qF  = uxF*uxF+uyF*uyF+uzF*uzF;
            const flow_float Un  = uxF*nx + uyF*ny + uzF*nz;
            const double TFd  = 0.5*(T0d+T1d);
            const double RFg  = thermo_R_mix(sp, ns, YF);
            const double cpFd = thermo_cp_mix(sp, ns, YF, TFd);
            const double cvFd = cpFd - RFg;
            const double gaFd = cpFd/cvFd;
            const flow_float c  = (flow_float)sqrt(gaFd*RFg*TFd);     // 凍結音速 (面平均 T, 面平均組成)
            // Ht・e は保存量から (EOS 厳密・datum 整合)
            const flow_float Ht = 0.5*( (roe_c[ic0]+Ps[ic0])/ro[ic0] + (roe_c[ic1]+Ps[ic1])/ro[ic1] );
            const flow_float eF = 0.5*( (flow_float)(roe_c[ic0]/ro[ic0]-0.5*q0d)
                                       +(flow_float)(roe_c[ic1]/ro[ic1]-0.5*q1d) );
            const flow_float r2E = 0.5*qF + eF - (flow_float)(cvFd*TFd);  // エントロピー波エネルギー成分

            // 接線基底
            flow_float ax, ay, az;
            if (fabs(nx) < 0.9) { ax=1.0; ay=0.0; az=0.0; } else { ax=0.0; ay=1.0; az=0.0; }
            flow_float t1x = ny*az - nz*ay, t1y = nz*ax - nx*az, t1z = nx*ay - ny*ax;
            const flow_float t1n = sqrt(t1x*t1x+t1y*t1y+t1z*t1z);
            t1x/=t1n; t1y/=t1n; t1z/=t1n;
            const flow_float t2x = ny*t1z - nz*t1y, t2y = nz*t1x - nx*t1z, t2z = nx*t1y - ny*t1x;

            // 固有値 (音響のみ c(') 込み) と TP スケーリング
            const flow_float cd  = (keepDissCprime == 1)
                                 ? lowMachCprime(c, sqrt(qF), Un, precondEps)
                                 : c;
            const flow_float lamA = fabs(Un) + cd;
            const flow_float lamU = fabs(Un);
            const flow_float sA = (flow_float)((double)roF/(2.0*gaFd*RFg));
            const flow_float sE = (flow_float)((double)roF/cpFd);
            const flow_float sS = (flow_float)((double)roF*TFd);

            // z_k = S_k |λ_k| (r_k · Δw) → D = Σ z_k r_k (エントロピー波のみ r_E=r2E)
            const flow_float ut1 = uxF*t1x + uyF*t1y + uzF*t1z;
            const flow_float ut2 = uxF*t2x + uyF*t2y + uzF*t2z;
            const flow_float rd1 = dw0 + (uxF-c*nx)*dw1 + (uyF-c*ny)*dw2 + (uzF-c*nz)*dw3 + (Ht-c*Un)*dw4;
            const flow_float rd2 = dw0 + uxF*dw1 + uyF*dw2 + uzF*dw3 + r2E*dw4;
            const flow_float rd3 = t1x*dw1 + t1y*dw2 + t1z*dw3 + ut1*dw4;
            const flow_float rd4 = t2x*dw1 + t2y*dw2 + t2z*dw3 + ut2*dw4;
            const flow_float rd5 = dw0 + (uxF+c*nx)*dw1 + (uyF+c*ny)*dw2 + (uzF+c*nz)*dw3 + (Ht+c*Un)*dw4;
            const flow_float z1 = sA*lamA*rd1;
            const flow_float z2 = sE*lamU*rd2;
            const flow_float z3 = sS*lamU*rd3;
            const flow_float z4 = sS*lamU*rd4;
            const flow_float z5 = sA*lamA*rd5;

            const flow_float coef = 0.5*sigmaF*sss;
            const flow_float Dro   = z1 + z2 + z5;
            const flow_float DroUx = z1*(uxF-c*nx) + z2*uxF + z3*t1x + z4*t2x + z5*(uxF+c*nx);
            const flow_float DroUy = z1*(uyF-c*ny) + z2*uyF + z3*t1y + z4*t2y + z5*(uyF+c*ny);
            const flow_float DroUz = z1*(uzF-c*nz) + z2*uzF + z3*t1z + z4*t2z + z5*(uzF+c*nz);
            const flow_float Droe  = z1*(Ht-c*Un) + z2*r2E + z3*ut1 + z4*ut2 + z5*(Ht+c*Un);

            res_ro_temp   -= coef*Dro;
            res_roUx_temp -= coef*DroUx;
            res_roUy_temp -= coef*DroUy;
            res_roUz_temp -= coef*DroUz;
            res_roe_temp  -= coef*Droe;
        }

        // スカラー輸送用の面質量流束 (連続式と整合: 散逸込みの総流束 ic0->ic1)
        // advGauge 時は参照流束 C∞·s を戻して物理流束にする (係数が場に依存する種輸送には
        // F∞ が telescoping しないためゲージ不可。種の free-stream 桁落ちは plan §8.2 の残課題)
        massflux[ip] = advGauge ? (res_ro_temp + CinfS) : res_ro_temp;

        atomicAdd(&res_ro[ic0]  , -res_ro_temp);
        atomicAdd(&res_roUx[ic0], -res_roUx_temp);
        atomicAdd(&res_roUy[ic0], -res_roUy_temp);
        atomicAdd(&res_roUz[ic0], -res_roUz_temp);
        atomicAdd(&res_roe[ic0] , -res_roe_temp);

        atomicAdd(&res_ro[ic1]  , res_ro_temp);
        atomicAdd(&res_roUx[ic1], res_roUx_temp);
        atomicAdd(&res_roUy[ic1], res_roUy_temp);
        atomicAdd(&res_roUz[ic1], res_roUz_temp);
        atomicAdd(&res_roe[ic1] , res_roe_temp);
    }

    __syncthreads();
}
