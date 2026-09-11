# HPC-01C：runtime 結束清理

日期：2026-09-09。基準工作樹 `d9eaea9`；2026-09-11 已由目前 HEAD 的
[F06 最終簽核](HPC_01C_F06_COMPLETION_REPORT.md)完成 HPC-01C。
證據根目錄：`outputs/hpc01/runtime-cleanup/`；最終彙整為 `acceptance.json`。

## 問題與結果

貼體 flow／transport runtime 原本只在 destructor 呼叫 PETSc destroy，忽略
回傳碼，呼叫端可能已列印完成摘要。新增 collective `Close()`：以固定順序
嘗試所有釋放，保存第一個返回錯誤，最後才共同決定成功／失敗。
重複 Close／destructor 不重新執行 destroy，失敗結果也不因重複 Close 而消失。
建構／求解已失敗時保留 noexcept 退棧清理。

CPU flow／VCA 在最後輸出關閉後明確 Close；graph／bifurcation 與 sequential
在最後 accepted history 建立後、CSV／completion manifest 發布前 Close。
因此清理失敗不列印成功摘要，耦合 CLI 亦不建立完成輸出目錄。
不改物理模型、PETSc 求解選項、時間步、partition 或任何既有檔案格式。

完整 lifecycle、物件順序與 embedding caller 責任見
[Runtime cleanup 契約](../architecture/RUNTIME_CLEANUP.md)。
本批的終止呼叫鏈稽核見 [入口稽核](HPC_01C_ENTRY_AUDIT.md)。

## 驗收

GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real64／Int32／MUMPS；
本機 1／2／3 ranks，OMP／BLAS 1 thread。所有相關 target 編譯完成，
無 compiler warnings，`git diff --check` 通過。

| 驗收 | 結果 |
|---|---|
| Runtime close 的逐物件返回錯誤 | flow／transport 各 12 個釋放位置，普通與空 rank fixture，world3 與 split1／2：144 次故障、144 次新 runtime 重試通過 |
| 清理完整性 | 每個 fault 都驗證所有釋放被嘗試一次、同群相同診斷、重複 Close 保留失敗、不重複釋放；retained PETSc references 回到單一測試 owner，之後完整釋放 |
| 建構失敗回歸 | 同次 `runtime_construction_failure_test` 的既有普通／空 rank／scalable options rollback 檢查全部返回 0 |
| 最終原生 CLI 範圍 | 112 作業：16 個前版參考、16 個健康、40 個故障、40 個新作業重試；168 份逐 rank 報告通過，無 timeout |
| CPU flow／VCA | 32 個 final velocity／pressure 比較 relative L2=0；另比對全部 74 份輸出（含 VTK／index／VCA manifest），逐位元一致 |
| Graph／sequential | 222 份健康／重試 CSV／manifest 與前版逐位元一致；28 個故障都沒有建立輸出目錄／completion manifest |
| 最終 graph subcommunicator | flow 的 explicit／Aitken 與 species graph，1+2 ranks 各比對 COMM_SELF；最大 relative L2 分別 `7.90598e-11`／`2.90249e-11`，通過既有 `1e-6` |

原生入口涵蓋 CPU flow／VCA、schema-v5 flow graph、schema-v6 species graph、
bifurcation executable、sequential explicit／fixed／Aitken；每個皆有 1／2 ranks。
在最後一個 rank，對 flow solver／Jacobian 或 transport solver／left matrix
的真實 destroy 返回後注入 PETSC_ERR_USER，逐 rank 核對 native status=1。
健康案例沿用原 nonlinear／mass gates；VCA 同時執行物種／reservoir 路徑。
正常及故障的 host peak RSS、wall time 與 solver stdout 保存在 rank reports。
本批是終止錯誤處理，不是效能工作，沒有宣稱組裝／求解加速或 RSS 改善。

## 證據版本與失敗嘗試

第一個 `native/` 在 graph 入口被拒絕：控制器先建立了輸出目錄，違反既有
「目錄必須不存在」規則。failed summary／logs 保留，修正控制器後重跑。

`native-final/` 的 112 個作業通過退出與數值比較，但中間版仍把耦合 Close
放在寫 manifest 之後。最終版將 Close 移到 writer 前，重建三個 coupling CLI
與 wrappers；`native-publication/` 再跑全部 80 個 coupling 作業並要求故障時
沒有輸出目錄。`acceptance.json` 保存 28 個中間版／最終版發布狀態對照。

最終原生證據選用 `native-final/` 中 source／binary 未再變動的 32 個 CPU
作業，加上 `native-publication/` 的 80 個 coupling 作業；不是把中間版的
coupling 結果當作最終版。最終 binary SHA-256 已逐一核對。
`cpu-output-comparison.json` 另補完整 CPU 檔案比較，沒有重跑未變更的求解。

## 重現

```bash
make -C solvers/cpu PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real \
  iga_navier_stokes runtime_cleanup_cli_test
make -C solvers/coupling petsc runtime_construction_failure_test \
  graph_cleanup_cli_test bifurcation_cleanup_cli_test sequential_cleanup_cli_test \
  multidomain_subcommunicator_test
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 timeout --kill-after=5s 240s \
  mpiexec --oversubscribe -np 3 solvers/coupling/runtime_construction_failure_test \
  /path/to/new-runtime-fixture
python3 scripts/hpc_runtime_cleanup_regression.py \
  --baseline-dir outputs/hpc01/runtime-cleanup/before --fixture-root outputs \
  --output-dir /path/to/new-native-output
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 timeout --kill-after=5s 180s \
  mpiexec --oversubscribe -np 3 solvers/coupling/multidomain_subcommunicator_test \
  /path/to/new-flow-groups
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 timeout --kill-after=5s 180s \
  mpiexec --oversubscribe -np 3 solvers/coupling/multidomain_subcommunicator_test \
  /path/to/new-species-groups species
```

CLI 控制器需要保留的 solver-stdout、checkpoint-write VCA、registry 及
sequential-output fixtures；case／binary hashes 與完整 argv 在 summaries。
預設完整執行；`--scope coupling` 僅供有明確既有 CPU 證據時重驗 coupling。

## 限制

返回錯誤注入發生在真實 PETSc destroy 完成後。這不是破壞 PETSc 內部狀態、
MPI 失聯或 rank 死亡測試，不能用來宣稱那些情況可原地恢復／無洩漏。
終止 Close 後不能再求解；重試建立新 runtime。GPU、完整 FSI 與前處理
未因本批重跑，因它們不使用這兩個 runtime 的新 Close 路徑。
F06 已由最終呼叫鏈與原生矩陣簽核；HPC-01D 亦已完成，HPC-03～09 的其餘項目
與跨節點驗收仍依各自任務追蹤。
