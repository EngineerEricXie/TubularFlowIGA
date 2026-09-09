# HPC-01C：記憶體報告關檔與 main 稽核

日期：2026-09-09。基準 `f9c7223` 加本批工作樹；HPC-01C 仍未完成。

## 問題與修正

`iga_solve --memory-report` 每次 Record 都 flush，最後由 ofstream destructor
關檔，未檢查關檔結果。使用既有 Linux
[text_close_preload](../../solvers/coupling/tests/text_close_preload.cpp)
讓指定報告的底層 fclose 在實際關閉後返回 EOF：修改前雙 rank CLI 的兩個程序
仍退出 0 並列印成功摘要，重現證據保留於 `reproduce/summary.json`。

[DistributedMemoryRecorder](../../solvers/cpu/include/MemoryReport.hpp) 新增 Close：
同群進入 `memory report close`，root 明確關閉並檢查結果，錯誤共同回報。
`iga_solve` 在 PETSc objects 清理後、最後成功摘要前呼叫它。未啟用報告時為
no-op；成功後重複 Close 無副作用；失敗後重複 Close 保留失敗且不再次關閉，
Record 在已關閉狀態下共同拒絕。Destructor 只作本地退棧後備。

Embedding caller 必須沿用相同 communicator／啟用設定並以相同順序呼叫
Record／Close。此 helper 不處理程序死亡或 MPI 內部失聯。關檔是晚期操作，
之前已寫出的場仍保留；本批不提供輸出回復或 checkpoint bundle 原子發布。

## 驗收

[原生控制器](../../scripts/hpc_memory_close_regression.py) 的 12 個作業、18 份
rank reports 全部符合預期。1／2 ranks 各包含修改前健康與故障重現、新版健康、
故障、重試、停用報告；新版故障的所有 ranks 返回 1，無 timeout，沒有最後成功
摘要。每份啟用報告都有原本六個 stages、完整 per-rank RSS arrays，格式未變。

10 個 transport 場比較 relative L2 均為 0，保留既有 `1e-6` 與零參考 `1e-12`
門檻；40 份 auxiliary outputs（fields、VTU、PVD、physiology JSON）逐位元相同。
故障案例也比較已產生的最終場，確認錯誤是報告關檔而非數值求解失敗。

[Memory report test](../../solvers/cpu/tests/test_memory_report.cpp) 另在三 rank world
及 split1+2 上通過 disabled／正常 Close、重複 Close、關檔錯誤重複回報、關閉後
Record 拒絕及新 recorder 重試。三個 communicator groups 各命中一次真實
fclose 回傳故障；同群 diagnostic 一致，分組各寫自己的報告。

環境：GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real64/int32；OMP／BLAS 1。
兩個受影響 binary 重建無新 compiler warning。以下單次健康作業資料只作
正確性驗收記錄，非效能比較；時間為 per-rank exclusive 最大值。

| Ranks | 組裝 s | 線性求解 s | 各 rank host peak bytes |
|---:|---:|---:|---|
| 1 | 0.002294 | 0.000842 | 39727104 |
| 2 | 0.002381 | 0.001743 | 39997440、39936000 |

## 重現與剩餘稽核

```bash
make -C solvers/cpu iga_solve memory_report_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
make -C solvers/coupling text_close_preload.so
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 timeout --kill-after=5s 60s \
  mpiexec --oversubscribe -np 3 env \
  LD_PRELOAD="$PWD/solvers/coupling/text_close_preload.so" \
  solvers/cpu/memory_report_test /path/to/new-unit-output split
python3 scripts/hpc_memory_close_regression.py \
  --baseline outputs/hpc01/memory-close/before/iga_solve \
  --fixture-root outputs/hpc01/transport-cli/cli-final \
  --output-dir /path/to/new-native-output
```

World unit 將最後一個參數換成 `world`，並用另一個新目錄。原生控制器需要先前
已接受的 transport CLI fixtures；數值輸入及 binary hashes 保存於 summary。
Ignored evidence 為 `outputs/hpc01/memory-close/`，包含 reproduction、build logs、
unit world／split logs、native summary／逐 rank 紀錄與 acceptance。

本輪也完成五個 main 的控制流程交界核對，見 [入口稽核](HPC_01C_ENTRY_AUDIT.md)。
其餘 coupling runner／runtime／adapter／executor 及共用 cleanup helper 的完整
內部鏈仍待 F06 簽核。此處沒有把任何 HPC-01C、HPC-03～09 任務勾選完成。
