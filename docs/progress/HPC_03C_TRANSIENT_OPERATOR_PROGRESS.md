# HPC-03C：完整固定幾何暫態 operator

日期：2026-09-09。基準 `b7059e0` 加本批修改。
狀態：完成 owned 固定幾何 backward-Euler operator；尚未完成暫態 Newton、
accepted clock、runtime commit／rollback 與 graph 接入，HPC-03C 保持未勾選。

## 實作

[ImmersedTransientDistributedOperator.hpp](../../solvers/cpu/include/ImmersedTransientDistributedOperator.hpp)
組合已驗證的分散式 history／凍結外力體積元件，以及共用 ghost、port、gauge、
owned Mat／Vec、required-state exchange 與全域診斷。沒有建立序列數值 runtime
來供給多 rank 的 block。

此 C++ 介面接受固定且 immutable 的 domain／volume／surface／ghost catalogs
及 geometry identity。這些 catalogs 在所有 trial time 都代表同一個 stationary
geometry，material wall velocity 定義為零。它不是接受任意移動 kinematics 後
把速度消除的 moving wrapper；既有 moving／evaluated-epoch 介面與限制不變。
舊序列 `ImmersedTransientFlowRuntime` 的 epoch-time guard 也未改動。

每個 owned cell 使用實際 dt、committed velocity history 與 frozen force
建立 conservative volume，加入 mixed trace 與既有 material-aware Nitsche
wall kernel。`wall_inertial_gamma0` 可為非零；ghost penalty、flow／pressure／
mean-normal-traction port 及 pressure gauge 沿用共用分散式積分與歸約。

[ImmersedStaticDistributedOperator.hpp](../../solvers/cpu/include/ImmersedStaticDistributedOperator.hpp)
新增本地 cell integrator 接點 `AssembleWithVolume`，原 `Assemble` 保留既有
steady integrator。共用的 Mat／Vec stashes、失敗協調、global port／wall-point／
pressure diagnostics 完整保留。`OwnedVolumeCells` 供 history／外力元件建立
同一份 work ownership；row ownership 與 integration ownership 仍可不同。

[ImmersedStaticFlowSetup.hpp](../../solvers/cpu/include/ImmersedStaticFlowSetup.hpp)
的新增 `FixedTransient` topology mode 保留 transient 的 numeric boundary-label
controller row 順序，diagnostic port vector 仍維持使用者宣告順序。steady 預設
保持原本 declaration-order rows。transient closed surface 採既有 transient
label completeness 政策，允許無 selected-wall point 的 positive cut cell；
steady 原有 closed-case 前置條件不變。新 mode 也納入 MPI configuration agreement。

共用 topology options 的 dt 仍為零；wrapper 只以新的 transient cell integrator
供給 volume／wall 係數，沒有將 steady solve 的 coupling clock 當作慣性項。
完整矩陣與殘差使用 `Freeze` 所發布的正 dt。

`Freeze` 先檢查完整 flow-target 相容性，再發布 frozen history／force／dt。
port values 可在 idle 依序供給，凍結期間拒絕修改；未知或 rank 間不同值仍由
共用協定拒絕。Release 取消 snapshot；Close 會嘗試關閉 inputs 與共用 operator
兩部分，即使第一部分回報錯誤也繼續釋放第二部分，最後回報第一個錯誤。

目前 `Diagnostics` 提供共用 physical assembly／port／wall-point／pressure
診斷；`ConservationDiagnostics` 是固定幾何 fluid flux／volume-divergence 視圖。
未宣稱已完整移植序列 transient 的所有 state hashes、wall-eta extrema 或
moving-continuity diagnostics，這些 runtime-level 整合仍待後續處理。

## 驗證

```bash
make -C solvers/cpu immersed_transient_distributed_operator_test \
  immersed_distributed_physics_test immersed_distributed_static_flow_test \
  CXX=mpicxx PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
python3 scripts/hpc_immersed_transient_operator_regression.py \
  --output-dir outputs/hpc03/transient-operator/final --split
python3 scripts/hpc_immersed_assembly_regression.py \
  --output-dir outputs/hpc03/transient-operator/static-flow \
  --kind physics --physics-mode flow --ranks 4
python3 scripts/hpc_immersed_assembly_regression.py \
  --output-dir outputs/hpc03/transient-operator/static-closed \
  --kind physics --physics-mode closed --ranks 2
python3 scripts/hpc_immersed_static_regression.py \
  --output-dir outputs/hpc03/transient-operator/static-newton --mode flow --ranks 2
make -C solvers/coupling iga_multidomain_flow CXX=mpicxx \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
```

本機 TsungYehLab、GCC 11.4、OpenMPI 4.1.2、PETSc 3.15.5 real64/int32；
OMP／BLAS threads=1。小型本機 correctness tests，MPI 使用允許 socket 的權限。
編譯開啟 warnings；受限編譯環境的 `opal_ifinit socket errno=1` 是環境訊息，
不是 compiler warning，授權 MPI runtime 已實際通過。

案例為 unit cube、4×1×1 cells、112 active nodes、compact quadrature。
flow、pressure、traction、closed、inertial（`wall_inertial_gamma0=0.6`）五模式
各跑 1／2／4 ranks，另跑 flow 的 split 1+2，共 16 MPI jobs、38 rank reports。
flow／inertial 為 451 rows，其餘為 449 rows。flow port declarations 故意反序，
驗證 controller row 仍與序列 `ImmersedActiveLayout` 一致。

每個 rank 的 COMM_SELF reference 僅供測試。非零 committed velocity、trial
velocity／pressure／controller／gauge 輸入共同激發矩陣與殘差；比較所有 owned
negative residual／Jacobian action 的 relative L2、全域 volume／surface／ghost
工作數、pressure measure、每個 port 的五項量測，以及 fluid conservation。
relative L2 門檻 1e-9，physical diagnostic 最大絕對差門檻 1e-12。

另檢查負值／rank 間不同的壁面慣性係數、未 Freeze／Close 後組裝與診斷拒絕、
frozen port 更新拒絕、incompatible all-flow targets 拒絕及補齊後重試。reference
的非零 wall temporal penalty 也明確檢查。這些是 operator input transactions，
不是 accepted timestep 的 runtime commit／rollback。

steady flow 4 ranks、closed 2 ranks 與 steady Newton flow 2 ranks 另外通過，
證明共享接點與預設 topology 沒有改變受測 steady 行為。正式 graph executable
也已重建；本批未重新執行全部 native graph，也未重跑未改變的長時間序列
closure suite。transient Newton 與完整 graph 模式仍未放行。

結果位於 `outputs/hpc03/transient-operator/`；`final/` 是最終 operator 驗收，
`initial/`／`full/` 是補最後配置與 Close guards 前的通過結果。資源報告包含
fixture／序列參考／啟動；operator 測試沒有 KSP solve，不能當成線性求解擴展
benchmark。加權分區沿用既有選項，本批 operator 案例使用預設 CellCount。

最終 operator 16 jobs／38 rank reports 全數 exit 0、無 timeout。最大 residual
relative L2 為 1.39734e-16，Jacobian-action relative L2
為 9.57066e-16，physical diagnostic 最大絕對差
1.77636e-15。最大 rank wall time 1.116142 s，
peak RSS 49020928 bytes。加上 steady 回歸，本批最終共
19 MPI jobs／46 rank reports 通過。原始 timings／RSS 保留於各 summary。
