# HPC-01C：Pressure-flow executor 的跨 rank 步驟邊界

日期：2026-09-08。HPC-01C 與整份清單仍未完成。
接續 [MPI 工具資產報告](HPC_01C_TOOL_ASSET_PROGRESS.md)。

## 問題與實作

輸入一致不代表執行過程中每個 rank 都會成功。原先
`PressureFlowComponentExecutor` 的殘差容器、port 結果、觀察 callback 或
replicated 0D runtime 若只在一個 rank 拋出例外，該 rank 會自行 abort，
其他 ranks 卻可能進入下一個 runtime 呼叫，甚至 finalize commit。

以修改前的 header、既有 FakeRuntime 與 3-rank MPI 重現：只讓 rank 1
的 precommit callback 失敗，得到 `failed_ranks=1 committed_upstream_ranks=2`。
這是部分提交的實際重現；FakeRuntime 本身不包含 MPI collective，因此
這個重現沒有宣稱已觀察到原生 PETSc deadlock。

新增不依賴 MPI 的 `PressureFlowExecutionSynchronization`，允許呼叫端
提供 operation outcome 與 all-converged callback。省略 policy 時維持
原本串行介面。MPI 呼叫端使用 `CollectivePressureFlowExecution`，借用
既有 communicator；沒有建立或釋放 communicator，也沒有使用固定 world。

每個可能失敗的本地資料階段會先協調結果，再前進；每個 runtime operation
完成或回報例外後，也會協調結果。每個 stage 最多呼叫一個 runtime operation，
其內部 collective 仍由 runtime 自己的錯誤協議負責。這層不能救回卡在
runtime collective 裡的 rank，也不能處理程序被作業系統殺死的情況。

主要邊界涵蓋初始 pressure／abort 記錄配置、BeginStep、輸入傳遞、SolveTrial、
port 讀取與驗證、edge／iteration 結果、relaxation、rollback、accepted ports、
precommit callback 與 PrepareCommitStep。Fixed／Aitken 只有在所有參與者
都同意時才收斂。FinalizeCommitStep 繼續遵守既有 noexcept 契約。

Abort 的 exception_ptr 陣列在 BeginStep 前配置；先逆序完成所有 runtime
的 abort，再建立文字診斷，避免報告字串配置失敗阻止後續清理。
保留 primary error，並附上各 domain 的 abort error。

原生 schema-v5 graph／bifurcation 與 sequential explicit 入口已接入 policy。
Executor 建構、precommit std::function 的配置及 accepted flow bookkeeping
加入 collective local stage。Sequential explicit 的 allocating pressure-map
參數也先完成協調，再呼叫 Advance。
`before_commit` 必須只執行本地工作；兩個已接入入口的 flow callback 符合此條件。
Graph 的 input／control 一致性檢查仍是同序執行的前提。

沒有更改殘差、Aitken 或數值求解公式，也沒有更改資料格式。
Accepted-result bookkeeping 失敗時會共同結束作業，但當時 runtime 已 commit；
本批沒有把它擴充為整個 committed step 的原地交易回復。

## 已取得的驗收證據

證據目錄：`outputs/hpc01/pressure-executor/`。

| 驗收 | 結果與範圍 |
|---|---|
| 舊版 defect | 3 ranks，rank 1 callback 失敗，其餘兩 ranks 仍提交 |
| 故障／健康重試 | 3-rank world，加上獨立 1／2-rank 子群；每組 20 個案例，共 60 次拒絕與 60 次成功重試 |
| 全域收斂 | 每組 Fixed、Aitken 各一次，只讓一個 rank 否決首次收斂；全部 ranks 共同 rollback 並在第二次 sweep 提交 |
| 清理再次失敗 | 同一次試驗中兩個 domain 的 abort 拋例外，仍完成全部 domain 的 abort，保留 primary 與兩份清理診斷 |
| 原生 3-rank graph | flow 與含 0D 的健康案例通過既有 Validate／ValidateZeroDClockRun；output 與 rank-local injection 設定差異共同退出 1 |
| 原生相容性 | MUMPS 設定下，兩個案例修改前後共 14 份 CSV 逐位元組相同 |
| Sequential MPI smoke | 原有 explicit、schema-v5 sequential、fixed、Aitken、subcycling 與故障案例通過；保留原門檻 |
| 完整 multidomain smoke | 延長 controller 上限後，原有全套案例退出 0；未改 fixture 或數值容許值 |
| 純 C++ 回歸 | `make -C solvers/coupling test` 通過，包含 pressure、species、graph、0D 與 FSI contracts |

20 類故障含實際 FakeRuntime solve／port／prepare／abort 與 precommit
callback 例外、rank-local missing pressure，以及在指定 operation outcome
注入 caught `std::bad_alloc`。後者驗證 outcome 協議，不宣稱測遍 C++
allocator 的每個呼叫點。重試使用同一 executor／runtime，在 abort 後重試
同一 step，並與串行 fixture 的結果精確比較。子群以不同次數額外執行 explicit
步驟，檢查誤用 world 的同步問題。

初次 native harness 漏帶既有 `hpc_failure_regression.py` 指定的 MUMPS 設定，
健康 runner 雖完成，既有 conservation gate 卻失敗並由測試 MPI_Abort(2)。
失敗紀錄保留在頂層 `native-healthy.log`；相同預設設定下的修改前 runner
得到完全相同的 7 份 CSV，見 `default-history-comparison.json`。
補回既定 PETSc 設定後的結果另存 `verified/`，沒有修改任何容許值。

完整 multidomain smoke 的第一次啟動在浸入式幾何計算期間達到 150 秒
controller timeout；到期前已確認該原生程序持續使用 CPU，退出碼 124 保留
於 `verified/acceptance-runs.json`，不算通過。確認原程序終止後，以相同
fixture／數值門檻、1800 秒上限重新執行，結果寫入 `smoke-final/`。
重跑已退出 0：sequential suite 約 22.5 秒，完整 multidomain suite 約
948.1 秒，包含原有浸入式、multi-island 與 schema-v6 species 驗收。
Species edge／domain／global 最大 residual 分別為 `6.7762635780344027e-20`、
`9.556334131140265e-15`、`8.890183375706226e-15`。這是 species executor
下一批修改之前的驗收；其來源與 binary 已封存，不能混同為下一批版本的證據。
`TUBULARFLOWIGA_KEEP_TEST_OUTPUT` 現在也適用於 multidomain smoke，供保留
實際案例與內部命令日誌；未設定時維持原本測試後移除案例的行為。

## 重現命令

環境：本機 WSL、GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real64／Int32；
OMP／BLAS 各 1。建置未新增 compiler warning。

```bash
make -C solvers/coupling pressure_flow_collective_failure_test \
  iga_multidomain_flow iga_1d_3d_explicit iga_1d_3d_bifurcation \
  multidomain_failure_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
make -C solvers/coupling test
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  timeout --kill-after=5s 90s mpiexec --map-by core --bind-to core -np 3 \
  solvers/coupling/pressure_flow_collective_failure_test
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps' \
  timeout --kill-after=5s 90s mpiexec --map-by core --bind-to core -np 3 \
  solvers/coupling/multidomain_failure_test /tmp/executor-healthy-new healthy
```

Native fixture 的目錄必須不存在。精確命令與退出碼另記於
`native-final/acceptance-runs.json`（6 筆最終 native 案例）與
`smoke-final/summary.json`；舊版 header／binary／重現程式保存在本批目錄。
測試 fixture 保留選項修改後，已重新建置 native 測試並重跑原生比較，
最終 14 份 CSV 的逐位元組比較見 `native-final/history-comparison.json`。

## 仍需補齊的邊界

下表為本批封存時的範圍；其後 species executor 的實作與獨立驗收見
[下一批報告](HPC_01C_SPECIES_EXECUTOR_PROGRESS.md)。

| 入口／元件 | 後續工作 |
|---|---|
| SpeciesPressureFlowComponentExecutor | Hydraulic iteration、donor／route 決策、transport accounting、precommit 與 abort 仍需等價協調；不可把本批 flow-only 覆蓋當成 species 覆蓋 |
| Graph runtime registry | Adapter／owner 集合與 registry 建構的配置／析構邊界仍需完整核對 |
| Sequential strong-fixed／Aitken CLI | 自行實作的強耦合迴圈沒有使用本批 executor，仍需檢查局部配置、history 與清理 |
| Sequential output | CSV／兩種 manifest 的 close 檢查與共同錯誤回報仍有缺口 |
| Flow runtime 與其他入口 | 接續既有報告做有限範圍的入口覆蓋稽核，再判定 HPC-01C 是否可勾選 |

本批沒有效能、GPU、跨節點或程序失聯恢復的驗收。新增加的同步成本尚未做
擴展性量測；單機小案例成功不能用來宣告 supercomputer scalability。

本批最終稽核使用 `source-pending.json`／`source-pending.tar.gz` 與
`pending-binaries/` 中已封存的版本；命名源自等待完整 smoke 期間。
後續 species 開發開始後，沒有把新的工作樹內容冒充本批受測來源。
