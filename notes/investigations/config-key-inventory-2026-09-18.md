# solverConfig キーの棚卸し (2026-09-18)

plan [config-key-pruning](../../plans/accepted/config-key-pruning.md) §5.1 #2 の分類表の素。
`solver_density_cuda/tools/config_key_inventory.py` の出力に、文書での言及先を足したもの。**判断は plan 側に書く**。

対象: 全 worktree の `case/**/solverConfig.yaml` **4063 本** (10 worktree, 内容が違う同名 0 本を別設定として計上)。
読み取り失敗 0 本。

列: **既定** = 実効既定値 (必須キーは既定なし)、**記載** = そのパスを書いている run 数、
**非既定** = 既定と違う値を書いている run 数 (既定が解決できないものは `—`)、**文書** = 言及があるファイル群。


## live な値キー (167)

| パス | 既定 | 記載 | 非既定 | 非既定の値 | 文書 |
| --- | --- | ---: | ---: | --- | --- |
| `bodyForce` | `?` | 165 | — |  | procedures/methods/plans(3) |
| `bodyForceCtrl` | `0` | 74 | 71 | 1 | procedures/methods/plans(1) |
| `bodyForceCtrlRelax` | `1.0` | 74 | 74 | 0.02 | procedures/methods/plans(1) |
| `bodyForceCtrlTarget` | `0.0` | 74 | 74 | 44.1248 | procedures/plans(1) |
| `detectNaN` | `0` | 186 | 186 | 1 | procedures/methods/plans(10) |
| `detectNaNInterval` | `1` | 186 | 176 | 10, 100, 50 | methods/plans(2) |
| `gpu` | `(required)` | 4060 | — |  | procedures/methods/plans(19) |
| `initial` | `(required)` | 4052 | — |  | procedures/methods/plans(5) |
| `keepDissCbCoeff` | `0.0` | 7 | 7 | 0.02, 0.1 | methods/plans(2) |
| `keepDissCbEps` | `0.10` | 7 | 0 |  | methods/plans(2) |
| `keepDissCoeff` | `0.05` | 204 | 94 | 0.015, 0.02, 0.03, 0.04 | methods/plans(4) |
| `keepDissCoeffMax` | `1.0` | 50 | 1 | 0.3 | methods/plans(1) |
| `keepDissCprime` | `1` | 12 | 7 | 0 | methods/plans(2) |
| `keepDissFdBlend` | `0` | 50 | 50 | 1 | methods/plans(2) |
| `keepDissJump` | `0` | 183 | 183 | 1, 2 | methods/plans(10) |
| `keepDissOpBlendRaw` | `0` | 41 | 41 | 1 | methods/plans(1) |
| `keepDissPrecond` | `0` | 167 | 162 | 1 | methods/plans(4) |
| `keepDissType` | `0` | 240 | 225 | 1, 2 | procedures/methods/plans(7) |
| `solver` | `(required)` | 4062 | — |  | procedures/methods/plans(41) |
| `condensation.condDTmaxStep` | `1.0` | 2 | 2 | 1.0e9 | procedures/methods/plans(1) |
| `condensation.condDgMaxStep` | `5.0e-3` | 2 | 2 | 1.0e9 | procedures/methods/plans(1) |
| `condensation.condEqDTmax` | `10.0` | 0 | 0 |  | methods/plans(2) |
| `condensation.condEqDgMax` | `0.05` | 0 | 0 |  | methods/plans(2) |
| `condensation.condEqRelax` | `1.0` | 0 | 0 |  | methods/plans(2) |
| `condensation.condEquilibrium` | `0` | 54 | 54 | 1, 2 | procedures/methods/plans(9) |
| `condensation.condEvapKelvin` | `0` | 0 | 0 |  | methods/plans(1) |
| `condensation.condEvapRmin` | `1.0e-9` | 0 | 0 |  | methods/plans(1) |
| `condensation.condEvaporation` | `1` | 14 | 11 | 0 | methods/plans(1) |
| `condensation.condFloat` | `1` | 345 | 11 | 0 | procedures/methods/plans(3) |
| `condensation.condGasSpecies` | `-1` | 743 | 743 | 0, 1 | procedures/methods/plans(5) |
| `condensation.condGrowthModel` | `0` | 690 | 151 | 1 | methods/plans(2) |
| `condensation.condGyarmathyC` | `3.18` | 4 | 4 | 0.0, 1.59, 12.72, 6.36 | methods/plans(2) |
| `condensation.condKantrowitz` | `0` | 703 | 692 | 1, 2, 3 | procedures/methods/plans(8) |
| `condensation.condKantrowitzGammaMode` | `0` | 3 | 3 | 1 | procedures/methods/plans(4) |
| `condensation.condLimiterMode` | `1` | 404 | 14 | 0 | procedures/methods/plans(5) |
| `condensation.condModel` | `0` | 758 | 744 | 1 | procedures/methods/plans(6) |
| `condensation.condN2LatentLowT` | `1` | 3 | 2 | 0 | procedures/methods/plans(1) |
| `condensation.condN2LiquidCp` | `2000.0` | 2 | 2 | 1500.0, 2500.0 | procedures/methods/plans(2) |
| `condensation.condN2PsatLowT` | `1` | 3 | 3 | 0 | procedures/methods/plans(1) |
| `condensation.condSigmaScale` | `1.0` | 95 | 95 | 0.97, 1.03, 3.0 | procedures/methods/plans(3) |
| `condensation.condSonicModel` | `-1` | 3 | 3 | 0 | procedures/methods/plans(4) |
| `condensation.condTwoTemp` | `0` | 2 | 2 | 1 | procedures/methods/plans(5) |
| `condensation.condVaporMassFraction` | `-1.0` | 33 | 33 | 0.7671 | procedures/methods/plans(2) |
| `condensation.condensation` | `0` | 873 | 707 | 1 | procedures/methods/plans(18) |
| `condensation.condensationSpecies` | `""` | 3 | 3 | H2O | procedures/methods/plans(2) |
| `condensation.nCondSpecies` | `0` | 870 | 870 | 1 | procedures/methods/plans(5) |
| `mesh.axisCentroidShift` | `1` | 1491 | 105 | 0 | methods/plans(8) |
| `mesh.axisRFloor` | `0.0` | 3 | 3 | 0.0003, 1e-12 | methods/plans(4) |
| `mesh.axisymMethod` | `0` | 4 | 4 | 1 | procedures/methods/plans(6) |
| `mesh.bndFirstOrder` | `0` | 102 | 92 | 1 | procedures/methods/plans(8) |
| `mesh.discretization` | `"cell"` | 2987 | 2603 | node | procedures/methods/plans(36) |
| `mesh.gradLSQ` | `0` | 80 | 38 | 1, 2 | procedures/methods/plans(7) |
| `mesh.hoopAreaFromClosure` | `0` | 20 | 11 | 1 | procedures/methods/plans(3) |
| `mesh.isAxisymmetric` | `0` | 2169 | 1499 | 1 | procedures/methods/plans(9) |
| `mesh.meshFileName` | `(required)` | 4062 | — |  | procedures/methods/plans(1) |
| `mesh.meshFormat` | `(required)` | 4062 | — |  | procedures/plans(1) |
| `mesh.nodeInletCornerWall` | `0` | 424 | 424 | 1 | procedures/methods/plans(4) |
| `mesh.nodeWallDirichlet` | `1` | 1412 | 4 | 0 | procedures/methods/plans(11) |
| `mesh.nodeWallStressEdgeKernel` | `1` | 3 | 2 | 0 | methods/plans(2) |
| `mesh.primPack` | `0` | 11 | 5 | 1 | procedures/methods/plans(2) |
| `mesh.renumber` | `"none"` | 14 | 14 | rcm | procedures/methods/plans(2) |
| `mesh.valueFileName` | `(required)` | 4037 | — |  | procedures/methods/plans(1) |
| `mesh.wallDistExtraPhysIDs` | `?` | 66 | — |  | procedures/plans(1) |
| `output.extraFields` | `?` | 139 | — |  | procedures/plans(2) |
| `output.level` | `1` | 748 | 686 | 0, 2 | procedures/methods/plans(8) |
| `physProp.Sc` | `0.7` | 79 | 2 | 1.0337e-05 | procedures/methods/plans(3) |
| `physProp.Sc_t` | `0.7` | 57 | 57 | 0.2, 0.9 | procedures/methods/plans(3) |
| `physProp.axisymMethod` | `0` | 1 | 1 | 1 | procedures/methods/plans(6) |
| `physProp.cp` | `(required)` | 4062 | — |  | procedures/methods/plans(22) |
| `physProp.gamma` | `(required)` | 4056 | — |  | procedures/methods/plans(26) |
| `physProp.isAxisymmetric` | `0` | 659 | 376 | 1 | procedures/methods/plans(9) |
| `physProp.isCompressible` | `(required)` | 4062 | — |  | methods/plans(1) |
| `physProp.pMin` | `1.0` | 113 | 112 | 1e-06, 20.0, 60.0 | procedures/methods/plans(6) |
| `physProp.prandtlLam` | `0.72` | 800 | 0 |  | procedures/plans(2) |
| `physProp.ro` | `(required)` | 4062 | — |  | procedures/methods/plans(43) |
| `physProp.roMin` | `1.0e-4` | 86 | 86 | 1e-06 | procedures/methods/plans(3) |
| `physProp.species` | `?` | 1544 | — |  | procedures/methods/plans(16) |
| `physProp.speciesDBFile` | `""` | 1268 | 1268 | species_db.yaml | procedures/methods/plans(4) |
| `physProp.speciesDiffusionMethod` | `1` | 238 | 58 | 0 | methods/plans(2) |
| `physProp.tMin` | `1.0e-4` | 86 | 86 | 1e-06 | procedures/methods/plans(3) |
| `physProp.thermCond` | `(required)` | 4062 | — |  | procedures/methods/plans(6) |
| `physProp.thermCondMethod` | `0` | 800 | 800 | 1 | procedures/plans(1) |
| `physProp.thermalMethod` | `(required)` | 4012 | — |  | procedures/methods/plans(15) |
| `physProp.thermoFloat` | `1` | 6 | 1 | 0 | procedures/methods/plans(3) |
| `physProp.thermoHrefTemp` | `0.0` | 1194 | 1191 | 298.15 | procedures/methods/plans(11) |
| `physProp.tracer` | `""` | 145 | 145 | exhaust | procedures/methods/plans(2) |
| `physProp.visc` | `(required)` | 4062 | — |  | methods/plans(7) |
| `physProp.viscMethod` | `(required)` | 4012 | — |  | procedures/methods/plans(10) |
| `physProp.chemistry.enabled` | `0` | 162 | 89 | 1 | procedures/methods/plans(2) |
| `physProp.chemistry.freezeBelowT` | `0.0` | 0 | 0 |  | procedures/methods/plans(1) |
| `physProp.chemistry.jacobianMode` | `1` | 162 | 155 | 2 | procedures/methods/plans(1) |
| `physProp.chemistry.mechanismFile` | `""` | 162 | 162 | mech.yaml, mech13.yaml | procedures/methods |
| `physProp.chemistry.tMaxReaction` | `6000.0` | 0 | 0 |  | procedures/methods/plans(1) |
| `space.convMethod` | `(required)` | 4062 | — |  | procedures/methods/plans(14) |
| `space.limiter` | `(required)` | 4037 | — |  | procedures/methods/plans(22) |
| `space.pRef` | `0.0` | 713 | 665 | 1000.0, 1000000.0, 101325.0, 1026.0 | procedures/methods/plans(8) |
| `space.roRef` | `0.0` | 292 | 286 | 1.1768, 1.177, 1.2 | procedures/methods/plans(1) |
| `space.uRef` | `?` | 292 | — |  | procedures/methods/plans(1) |
| `time.bdfOrder` | `2` | 383 | 50 | 1 | methods/plans(3) |
| `time.dualTime` | `(required)` | 4004 | — |  | procedures/methods/plans(7) |
| `time.nStepInner` | `1` | 3883 | 3427 | 10, 100, 1000, 10000 | procedures/methods/plans(16) |
| `time.nSubIterDualTime` | `20` | 473 | 266 | 10, 12, 13, 15 | procedures/methods/plans(5) |
| `time.outStepInterval` | `(required)` | 4062 | — |  | procedures/methods/plans(5) |
| `time.outStepStart` | `(required)` | 4008 | — |  | — |
| `time.timeIntegration` | `(required)` | 4054 | — |  | procedures/methods/plans(18) |
| `time.unsteady` | `(required)` | 4004 | — |  | procedures/methods/plans(12) |
| `time.deltaT.axisTimestepBeta` | `0.0` | 39 | 39 | 2.0 | procedures/methods/plans(3) |
| `time.deltaT.blockDPLUR` | `0` | 3616 | 3328 | 1 | procedures/methods/plans(12) |
| `time.deltaT.blockDPLURDiagCache` | `0` | 2 | 2 | 1 | procedures/plans(2) |
| `time.deltaT.blockDPLURDqPack` | `0` | 9 | 5 | 1 | procedures/methods/plans(2) |
| `time.deltaT.cfl` | `(required)` | 4062 | — |  | procedures/methods/plans(40) |
| `time.deltaT.cfl_pseudo` | `(required)` | 4004 | — |  | procedures/methods/plans(34) |
| `time.deltaT.condRealizProject` | `1` | 2 | 2 | 0 | methods/plans(2) |
| `time.deltaT.control` | `(required)` | 4062 | — |  | procedures/methods/plans(6) |
| `time.deltaT.detectNaN` | `0` | 2885 | 2787 | 1 | procedures/methods/plans(10) |
| `time.deltaT.detectNaNInterval` | `1` | 0 | 0 |  | methods/plans(2) |
| `time.deltaT.dt` | `(required)` | 4062 | — |  | procedures/methods/plans(32) |
| `time.deltaT.dt_max` | `(required)` | 4052 | — |  | procedures |
| `time.deltaT.dt_min` | `(required)` | 4052 | — |  | procedures/plans(1) |
| `time.deltaT.ducrosLimiter` | `0` | 56 | 1 | 1 | procedures/plans(2) |
| `time.deltaT.implicitRelax` | `1.0` | 1859 | 1013 | 0.1, 0.3, 0.4, 0.5 | procedures/methods/plans(14) |
| `time.deltaT.implicitSolvePrecision` | `0` | 30 | 5 | 1 | methods/plans(4) |
| `time.deltaT.lineDtDirectional` | `0` | 35 | 35 | 1 | procedures/methods/plans(1) |
| `time.deltaT.lineImplicit` | `0` | 71 | 71 | 1 | procedures/methods/plans(4) |
| `time.deltaT.lineKFreeze` | `0` | 35 | 34 | 1 | procedures/methods/plans(4) |
| `time.deltaT.lineViscCoupling` | `0` | 34 | 34 | 1 | procedures/methods/plans(1) |
| `time.deltaT.lineViscousDtRelief` | `0.0` | 33 | 33 | 1.0 | procedures/methods/plans(3) |
| `time.deltaT.lowMachPrecond` | `0` | 2192 | 199 | 1, 2, 3 | procedures/methods/plans(14) |
| `time.deltaT.lowMachThornber` | `0` | 7 | 4 | 1 | procedures/methods/plans(2) |
| `time.deltaT.monitorInterval` | `1` | 433 | 433 | 10, 100, 2, 5 | procedures/methods/plans(5) |
| `time.deltaT.multispeciesRhoYCommonLimiter` | `0` | 25 | 4 | 1 | procedures/methods/plans(2) |
| `time.deltaT.passiveFct` | `1` | 131 | 34 | 0 | procedures/methods/plans(2) |
| `time.deltaT.passiveFctPrelimit` | `0` | 0 | 0 |  | procedures/plans(1) |
| `time.deltaT.passiveFctSweeps` | `100` | 0 | 0 |  | procedures/plans(2) |
| `time.deltaT.passiveFctTol` | `1.0e-6` | 2 | 2 | 1e-30 | procedures/plans(2) |
| `time.deltaT.passiveFctTolAbs` | `1.0e-30` | 2 | 2 | 0.0, 1e+30 | procedures/plans(2) |
| `time.deltaT.passiveImplicitCoupling` | `-1` | 142 | 142 | 0, 1 | procedures/methods/plans(1) |
| `time.deltaT.passiveImplicitRelax` | `-1.0` | 11 | 11 | 1.0 | procedures/methods/plans(1) |
| `time.deltaT.passiveScalarScheme` | `1` | 440 | 14 | 0 | procedures/methods/plans(2) |
| `time.deltaT.precondEps` | `0.15` | 184 | 24 | 0.03, 0.05, 0.08, 0.1 | procedures/methods/plans(5) |
| `time.deltaT.scalarCflMax` | `-1.0` | 1 | 1 | 2.0 | procedures/methods/plans(1) |
| `time.deltaT.speciesFaceReconstruction` | `0` | 431 | 288 | 1, 2 | procedures/methods/plans(5) |
| `time.deltaT.speciesImplicitCoupling` | `0` | 498 | 461 | 1, 2 | procedures/methods/plans(6) |
| `time.deltaT.speciesImplicitRelax` | `1.0` | 10 | 0 |  | procedures/methods/plans(1) |
| `time.deltaT.updateGuardAlpha` | `0.0` | 11 | 10 | 0.5 | procedures/methods/plans(3) |
| `time.last.control` | `(required)` | 4062 | — |  | procedures/methods/plans(6) |
| `time.last.nStepOuter` | `(required)` | 4062 | — |  | procedures/methods/plans(4) |
| `turbulence.dilatationCorrection` | `2` | 1977 | 409 | 0, 1 | procedures/methods/plans(6) |
| `turbulence.kInit` | `0.0` | 1 | 1 | 1.0 | procedures/methods/plans(3) |
| `turbulence.katoLaunder` | `0` | 965 | 939 | 1 | procedures/methods/plans(7) |
| `turbulence.model` | `(required)` | 2775 | — |  | procedures/methods/plans(17) |
| `turbulence.nodeKwfDirichlet` | `1` | 3 | 2 | 0 | methods/plans(3) |
| `turbulence.nodeOmegaWfDirichlet` | `0` | 346 | 314 | 1 | methods/plans(5) |
| `turbulence.omegaInit` | `0.0` | 1 | 1 | 1000.0 | procedures/methods/plans(3) |
| `turbulence.scalarDiffusion` | `1` | 1984 | 8 | 0 | procedures/methods/plans(3) |
| `turbulence.sstCrossDiffJac` | `0` | 72 | 72 | 1 | procedures/plans(2) |
| `turbulence.sstEnergyIncludesK` | `0` | 38 | 35 | 1 | procedures/methods/plans(6) |
| `turbulence.sstEnergyKSource` | `0` | 3 | 3 | 1 | procedures/methods/plans(3) |
| `turbulence.sstEnergyWallFunction` | `0` | 5 | 5 | 1 | methods/plans(3) |
| `turbulence.sstIsotropicStress` | `0` | 3 | 3 | 1 | procedures/methods/plans(3) |
| `turbulence.sstNodeWallKPin` | `1` | 2 | 2 | 0 | procedures/methods/plans(2) |
| `turbulence.sstOmegaProdFromPk` | `1` | 8 | 5 | 0 | procedures/methods/plans(4) |
| `turbulence.sstSigmaBlend` | `1` | 8 | 5 | 0 | procedures/methods/plans(4) |
| `turbulence.sstThermalWallFunction` | `0` | 19 | 18 | 1, 3 | procedures/methods/plans(4) |
| `turbulence.turbulentPrandtl` | `0.85` | 578 | 578 | 0.9 | procedures/methods/plans(2) |
| `turbulence.turbulentSchmidt` | `0.7` | 163 | 39 | 0.2, 0.5 | procedures/methods/plans(1) |
| `turbulence.wallTreatmentSST` | `1` | 1385 | 638 | 0 | procedures/methods/plans(21) |

## 拒否専用の参照 (12) — 書くと起動時に落ちる

| パス | 既定 | 記載 | 非既定 | 非既定の値 | 文書 |
| --- | --- | ---: | ---: | --- | --- |
| `mesh.nodeAxisDirichlet` | `None` | 414 | — |  | procedures/methods/plans(11) |
| `mesh.nodeAxisUrDirichlet` | `None` | 21 | — |  | procedures/methods/plans(2) |
| `mesh.nodeMidpointFx` | `None` | 36 | — |  | procedures/methods/plans(2) |
| `mesh.nodeReconEdgeMidpoint` | `None` | 24 | — |  | procedures/methods/plans(2) |
| `mesh.nodeValueAtNode` | `None` | 73 | — |  | procedures/methods/plans(5) |
| `time.dualTime_InnerLoop` | `None` | 0 | — |  | plans(1) |
| `time.nInnerLoop` | `None` | 0 | — |  | procedures/plans(1) |
| `time.last.nStep` | `None` | 0 | — |  | procedures/plans(1) |
| `turbulence.DESmode` | `None` | 63 | — |  | methods/plans(3) |
| `turbulence.LESmodel` | `None` | 1259 | — |  | methods/plans(4) |
| `turbulence.LESorRANS` | `None` | 1259 | — |  | procedures/methods/plans(13) |
| `turbulence.RANSmodel` | `None` | 904 | — |  | procedures/methods/plans(4) |

## 起動時拒否テーブル removed[] (8)

| パス | 既定 | 記載 | 非既定 | 非既定の値 | 文書 |
| --- | --- | ---: | ---: | --- | --- |
| `mesh.gradLSQDegenThresh` | `None` | 0 | — |  | procedures/methods/plans(2) |
| `time.deltaT.implicitRelaxSST` | `None` | 4 | — |  | procedures/plans(4) |
| `time.deltaT.lineDtWallRelief` | `None` | 1 | — |  | procedures/plans(2) |
| `turbulence.C_DES_ke` | `None` | 0 | — |  | procedures/methods/plans(2) |
| `turbulence.C_DES_kw` | `None` | 0 | — |  | procedures/methods/plans(2) |
| `turbulence.wmlesNewtonMaxIt` | `None` | 0 | — |  | procedures/methods/plans(2) |
| `turbulence.wmlesNewtonTol` | `None` | 0 | — |  | procedures/methods/plans(2) |
| `turbulence.wmlesPrt` | `None` | 0 | — |  | procedures/methods/plans(2) |

## 節そのもの (10) — キーではない

| パス | 既定 | 記載 | 非既定 | 非既定の値 | 文書 |
| --- | --- | ---: | ---: | --- | --- |
| `condensation` | `None` | 0 | — |  | procedures/methods/plans(18) |
| `mesh` | `None` | 0 | — |  | procedures/methods/plans(45) |
| `output` | `None` | 0 | — |  | procedures/methods/plans(11) |
| `physProp` | `None` | 0 | — |  | procedures/methods/plans(11) |
| `space` | `None` | 0 | — |  | procedures/methods/plans(6) |
| `time` | `None` | 0 | — |  | procedures/methods/plans(33) |
| `turbulence` | `None` | 0 | — |  | procedures/methods/plans(41) |
| `physProp.chemistry` | `None` | 0 | — |  | procedures/methods/plans(10) |
| `time.deltaT` | `None` | 0 | — |  | procedures/methods/plans(8) |
| `time.last` | `None` | 0 | — |  | procedures/plans(2) |
