# HPC-03C：固定幾何暫態 Newton 與交易時鐘

日期：2026-09-09。基準 `eb69117` 加本批修改。
狀態：分散式 C++ runtime 已通過驗收；固定幾何暫態 case／graph 尚待接入，
HPC-03C 保持未勾選。

## 實作與範圍

[ImmersedTransientDistributedRuntime.hpp](../../solvers/cpu/include/ImmersedTransientDistributedRuntime.hpp)
將已完成的 owned transient operator 接上共用 MPI Newton／KSP。owned committed、
prepared、candidate vectors 與 velocity history 均保持分散；沒有在 production
runtime 複製完整序列場。geometry catalogs 仍為 replicated immutable metadata。

[ImmersedDistributedNewtonRuntime.hpp](../../solvers/cpu/include/ImmersedDistributedNewtonRuntime.hpp)
抽出既有 steady Newton 引擎，保留預設 GMRES、右預條件 LU／MUMPS、真實線性
殘差、global block norms、flow-controller gates 與 backtracking。steady prefix
仍為 `immersed_static_`，新 runtime 使用 `immersed_transient_`；prefix 與 PETSc
options 在 communicator 內必須一致。舊 diagnostic 型別與錯誤階段名稱保留。

`BeginTrial` 檢查所有 rank 的 source／target 時間、步數與正 dt，凍結上一步速度
及外力成功後才發布 trial clock。Newton 失敗恢復 committed field，保留凍結輸入
供相同 trial 重試。`Rollback` 同樣保留 history；`AbortTrial` 才會丟棄 history
並回到 accepted clock。外部 coupling adapter 重做邊界條件時必須使用後者。

`PrepareCommit` 只準備 owned candidate；`FinalizeCommit` 才同時發布場與時鐘，
並釋放 frozen history。取消 prepare、abort、重複 finalize 或 Close 後 finalize
均不推進時鐘。`AbortTrial` 不恢復 caller 設定的 port controls，graph adapter
仍須保存與恢復 accepted controls。所有 MPI 方法為 collective；noexcept
publication 方法由所有成員在共同決策後呼叫。

`Diagnostics()` 描述 Newton phase，`Clock()` 描述更長的 frozen-trial lifetime；
Begin 後及失敗重試期間兩者的 `trial_active` 可不同。守恆視圖沿用固定幾何 fluid
flux／volume divergence。未宣稱完整移植序列 transient 的 hashes、wall-eta
extrema 或 moving-geometry diagnostics；既有檔案格式與 epoch guards 保留。

## 抵消造成的序列守恆誤判

第一個 flow regression 在序列參考的守恆診斷失敗。進出口約 ±1e-4 的通量先在
同一 aggregate 內抵消，舊 algebraic identity gate 只剩近零的 net-flow scale，
因此將正常加總 roundoff 誤判為不一致。

`ImmersedTransientFlowRuntime` 現在保留 surface fluid 與 material-wall 原始
quadrature 通量的絕對值和、累加次數，使用 gamma-n 累加誤差界與既有 epsilon
allowance 檢查 identity。moving owner 的二次檢查使用同一份原始統計。
物理 continuity、divergence、leakage 公式與 normalization 沒有改動。

聚焦測試以三個 ±1e-4／1e-20 contributions 重現舊誤判，實測差異
−3.55253e-21、roundoff allowance 6.21725e-18；兩倍 allowance 的不一致仍被拒絕。
另驗證舊 material cancellation、非法／非有限統計、累加界失效，以及 legacy
moving owner 的 stationary solve、inner／outer 有限缺陷注入與 commit 保留紀錄。
此聚焦測試首次提供錯誤的 target epoch 區間而被既有 guard 拒絕；修正 fixture
為 `[0, 0.125]` 後通過，production epoch guard 未改動。

## 驗證命令與門檻

```bash
make -C solvers/cpu immersed_transient_distributed_runtime_test \
  immersed_conservation_roundoff_test immersed_distributed_static_flow_test \
  CXX=mpicxx PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
python3 scripts/hpc_immersed_transient_runtime_regression.py \
  --output-dir outputs/hpc03/transient-runtime/final --split
python3 scripts/hpc_immersed_transient_runtime_regression.py \
  --output-dir outputs/hpc03/transient-runtime/weighted \
  --ranks 4 --modes flow empty-work --weighted
python3 scripts/hpc_immersed_static_regression.py \
  --output-dir outputs/hpc03/transient-runtime/static-initial --ranks 2 --mode flow
python3 scripts/hpc_immersed_static_regression.py \
  --output-dir outputs/hpc03/transient-runtime/static-closed --ranks 2 --mode closed
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  mpiexec -np 1 solvers/cpu/immersed_conservation_roundoff_test
make -C solvers/coupling iga_multidomain_flow CXX=mpicxx \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
```

TsungYehLab、GCC 11.4、OpenMPI 4.1.2、PETSc 3.15.5 real64/int32，OMP／BLAS
threads=1。MPI 以允許 socket 的本機權限執行；無跨節點或 CUDA 驗收。編譯開啟
warnings，受限編譯中的 `opal_ifinit socket errno=1` 為環境訊息。

unit cube 的 4×1×1 cells／112 nodes，empty-work 使用 2×1×1；compact
quadrature。flow、pressure、closed、inertial（gamma=0.6）、empty-work 各測
1／2／4 ranks，另測 flow split 1+2，以及 weighted flow／empty-work 4 ranks。
每個案例接受兩步 dt=0.125，第二步改變 flow targets；closed 案例使用非零初值
與外力。序列 reference 只存在測試內，按固定 global node IDs 重綁 target epoch。

預先設定 velocity／pressure／scalar 各區塊 relative L2 < 1e-6；reference norm
< 1e-10 的近零區塊改驗 absolute L2 < 1e-10。port flow／pressure／traction、
volume divergence／surface／wall flux 的比較門檻為 1e-11 + 1e-6 |reference|。
每次 linear solve 要求正收斂原因、true relative residual < 1e-6。第二步另以
清零 previous velocity 的 reference 求解，要求 velocity 最大差 > 1e-9，
防止測試只驗到重複單步求解。

負面案例包括單 rank 錯誤時鐘、rank 間不同但各自合法的 dt、凍結外力失敗、
凍結期間修改輸入、線性更新後的候選失敗、prepare failure、取消 prepare、
abort／re-Begin、重複 finalize 及 Close 後禁止發布。重試不得重新評估 frozen force。

## 結果與下一步

最終 18 個 transient MPI jobs／46 rank reports 全部通過，無 timeout。
最大 field relative L2 1.47572e-9；physical comparison 使用容許值的最大比例
6.3501e-8。history discriminator 最小差 1.24364e-7，明確超過 1e-9 門檻。
最大 rank wall time 7.290032 s、peak RSS 62,324,736 bytes；每 rank 累計
assembly 最大 3.41442 s，linear solve 最大 0.12358 s。

steady flow／closed 各 2 ranks 與聚焦 moving-owner 檢查亦通過，總共
21 個有效 MPI jobs／51 rank reports。正式 native graph executable 重建成功，
無新增 compiler warning；本批沒有重新執行完整 native graph 矩陣。

timing 包含各次 retry 的 runtime 工作；wall／RSS 另含 fixture、序列參考、
MPI 啟動。部分 correctness jobs 與編譯／其他小測試並行，未做固定核心的
獨立效能重複，因此不能宣稱 scaling 或加權分區加速。

證據位於 `outputs/hpc03/transient-runtime/`，包含 controller summaries、每 rank
timing／RSS／日誌、build logs 與 archive binaries。`initial` 保存原 roundoff
誤判；`conservation-final` 保存 epoch fixture 錯誤；有效聚焦結果位於
`conservation-epoch-fixed`。本批未重跑長時間完整 moving closure suite。

下一步為 transient graph adapter、case 配置與正式多步入口：coupling retry
需從同一 accepted velocity 重新凍結，提交前要同步所有 domain 的場與時鐘。
只有正式 graph 的 MPI 比較、守恆、rollback／precommit failure 通過後，才對
已支援固定幾何模式移除 size-1 限制並完成 HPC-03C。
