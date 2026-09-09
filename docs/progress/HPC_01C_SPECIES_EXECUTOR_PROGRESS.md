# HPC-01C：Species executor 的錯誤與排程協調

日期：2026-09-08。HPC-01C 及整份清單仍未完成。
接續 [Pressure-flow executor 報告](HPC_01C_PRESSURE_EXECUTOR_PROGRESS.md)。

## 問題與修正

`SpeciesPressureFlowComponentExecutor` 原先在一個 rank 的本地配置、
transport accounting 或 precommit callback 失敗時，自行執行 abort，
其他 ranks 可能繼續求解或提交。若 ranks 算出不同 species donor，還可能
進入不同的 transport domain 順序。

使用修改前的 header 與既有 FakeStagedRuntime，在 3 ranks 重現：
rank 1 的 precommit callback 失敗，得到 `failed_ranks=1 committed_ranks=2`。
另一個案例在 rank 1 反轉相容的兩端流量；舊版所有 ranks 都接受，但第一個
transport domain 不同，`schedule_min=0 schedule_max=1`。
這兩個重現使用本地 fake runtime，證明部分提交與排程分歧，並不宣稱觀察到
原生 PETSc collective 卡死。

新增可選、與 MPI 無關的 `SpeciesPressureFlowExecutionSynchronization`。
MPI 呼叫端經 `CollectiveSpeciesPressureFlowExecution` 借用 communicator，
提供 outcome、all-converged 與 same-schedule 三個 callback；必須全部提供
或全部省略。省略時維持既有串行使用方式。

已加入以下邊界：

- 初始 pressure 與 abort 記錄配置、BeginStep、hydraulic trial／rollback、
  port 讀取與驗證、edge／iteration 結果、relaxation 與全群收斂。
- Donor 選擇、route 建立、拓樸排序與 schedule 序列化。以長度框定的字串
  比對供應端 key／方向及 transport domain 順序；分歧在任何 transport
  solve 前共同拒絕。合法且一致的反向流仍會使用反向 donor。
- Donor concentration 讀取、本地輸入集合、SetTransportConcentration、
  transport solve、native accounting 讀取、本地 catalog 及所有 amount gates。
- 最終 transport port 觀察、precommit callback、PrepareCommitStep，以及
  逆序 abort。錯誤記錄先配置；所有 abort 完成後才組合文字，保留 primary
  與各 domain 的清理錯誤。

新增 `SpeciesPressureFlowStepResult::transport_ports` 保存最終完整 port
狀態。原生 schema-v6 callback 改為讀取已完成的本地結果，不再在 callback
中逐一呼叫可能含 collective 的 runtime getter。原生 accepted-result
bookkeeping 也加入 collective local stage，並將 transport map 移入既有
輸出記錄，維持原來的五種 CSV 格式。
Donor ownership 仍僅在全部 Prepare／Finalize 成功後提交；失敗 trial 不會
更新已提交的 donor 記錄。

每個 stage 最多呼叫一個 runtime operation。Runtime 必須處理自身內部的
collective 錯誤；外層只能協調已返回的 outcome，不能救回卡在內部 collective
的 rank。`before_commit` observer 必須只執行本地工作。Graph、controls 與
species metadata 的既有啟動一致性檢查仍是前提。
沒有改動 hydraulic／species 公式、物種單位、tolerance 或檔案格式。

## 驗收

本批 ignored 證據：`outputs/hpc01/species-executor/`。

| 驗收 | 結果與範圍 |
|---|---|
| 舊版重現 | 3 ranks 的部分提交與不同 transport 順序皆重現 |
| 新 MPI 故障／重試 | world 3 ranks：36 組；獨立 1-rank 子群：35 組；2-rank 子群：36 組；合計 107 次拒絕、107 次健康重試 |
| 全群收斂 | 三種 communicator 配置各測 Fixed／Aitken；單一 rank 否決首輪時，全部共同 rollback 並在第二輪提交 |
| 合法反向流 | 三種配置合計 6 次；驗證反向 domain 順序及兩種 species 的接收濃度 |
| 排程分歧 | 2／3-rank 案例在 schedule agreement 被拒絕，transport solve 次數為 0，沒有提交 donor ownership |
| 原生 schema-v6 smoke | 執行原有 species 子集：1／2 ranks、graph permutation、precommit failure、非法 binding、非法 1D mode、非 prescribed velocity 及既有 amount gates |
| 新舊相容性 | 1／2-rank 案例各 5 份 CSV，合計 10 份與修改前逐位元組相同 |
| 純 C++ contracts | `make -C solvers/coupling test` 通過，含 pressure、species、graph、0D 與 FSI contracts |

故障包含實際 fake runtime 的 transport／prepare／abort／port／accounting
例外、missing pressure、標準與非標準 callback 例外，以及在指定 outcome
注入 caught `std::bad_alloc`。最後一類驗證協調協議，不代表逐個測遍 allocator
呼叫點。每次失敗確認所有 domain 未 commit、donor ownership 未提交、各 ranks
取得相同診斷；再於相同 executor／runtime 重試同一 step，與串行 fixture
精確比對 port、routing 及 edge／domain／global amount 結果。

Native 最大 edge／domain／global residual 分別為
`6.7762635780344027e-20`、`9.556334131140265e-15`、
`8.890183375706226e-15`，與上一批受測版本相同。
MPI 與 native controllers 有 timeout，預期故障須有限時間返回，timeout
不算通過。本批沒有執行新的全套浸入式 smoke；上一批全套 smoke 的完成版本
與 948.1 秒紀錄保留在 pressure-executor 目錄，不冒充本批受測來源。

初次 MPI 測試 build 因將 serial `main` 改名作為 fixture 而出現 return-type
warning；已在原 serial test 加上明確 `return 0`，最終 build 沒有新增 compiler
warning。初始 build log 保留，沒有把原始警告藏成成功紀錄。

## 重現命令

環境：本機 WSL、GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real64／Int32；
OMP／BLAS 各 1。保持既有的每案例 PETSc 求解器設定。

```bash
make -C solvers/coupling petsc species_collective_failure_test \
  multidomain_flow_smoke_test test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  timeout --kill-after=5s 120s mpiexec --map-by core --bind-to core -np 3 \
  solvers/coupling/species_collective_failure_test
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  TUBULARFLOWIGA_KEEP_TEST_OUTPUT=1 \
  timeout --kill-after=5s 180s solvers/coupling/multidomain_flow_smoke_test --species-only
```

`--species-only` 只選取既有 smoke 的 species 部分，未更改其中的 fixture 或
驗收標準。不帶此選項時仍執行原本全套案例；也可用新增的
`make -C solvers/coupling multidomain-species-test` 執行子集。

## 接續工作

HPC-01C 接續核對 graph registry／adapter owner 集合的配置與析構、sequential
strong-fixed／Aitken 迴圈、sequential CSV／manifest close 與剩餘 flow runtime
入口。完整 checkpoint 發布與跨 rank 數續跑仍屬 HPC-05。
Accepted-result bookkeeping 失敗可共同終止作業，但 runtime 此時可能已 commit；
本批沒有新增 committed step 的原地回復。也沒有程序失聯、GPU、跨節點或
擴展性驗收，新同步成本尚未量測。清單保持 9／38，不宣告整份 goal 完成。

最終執行記錄、執行當下 binary 雜湊及 CSV 比較位於 `accepted/`；
`source-final.json`／`source-final.tar.gz` 與 `after-binaries/` 保存最終受測版本，
`acceptance.json` 記錄稽核結果。`before/` 保存舊 header、runner 與來源快照。
