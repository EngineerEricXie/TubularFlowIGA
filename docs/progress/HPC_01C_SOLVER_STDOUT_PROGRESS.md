# HPC-01C：CPU solver 的 stdout 錯誤邊界

日期：2026-09-08。基準 commit：`df535cf`；HPC-01C 仍未完成。
本批 CPU solver stdout 邊界實作與驗收完成；其餘 HPC-01C 邊界繼續追蹤。
證據：`outputs/hpc01/solver-stdout/`。

## 修改

CPU flow、configured transport 與 legacy transport 的輸入摘要、restart／checkpoint、
VTKHDF 初始化與完成摘要，在原有共同本地階段內明確 flush 並檢查 fail/bad 狀態。
Flow runtime 的元素組裝摘要、Newton iteration／convergence、outlet iteration／state
亦同。只在主執行緒輸出，不改 OpenMP worker、數值模型或 PETSc 插入順序。

Flow／configured transport 的最後 phase profile 原本位於主 try/catch 外；現在
使用共同本地階段，輸出失敗使全部 rank 退出 1，然後共同 finalize PETSc。
串流例外 mask 與 ownership 不變，沒有把 default ostream 改成全域 throwing 模式。

`test_tool_stdout_failure.cpp` 延伸為三個原生 CLI wrapper，攔截 rank 0 的真實
stdout buffer；依指定訊息片段注入 short-write、flush、runtime_error、bad_alloc
或非標準例外。控制器逐 rank 核對 native exit、診斷階段、timeout 與 RSS report，
每次故障之後以新作業重試，比較 final field；flow 同時比較 velocity 與 pressure。

新增 `_test` 執行檔與 FSI exporter 的 ignore 規則，避免測試產物混入 commit。
不改或移除 test sources。

## 案例與限制

本機 Open MPI 4.1.2、PETSc 3.15.5 real64／Int32／MUMPS，GCC 11.4。
1／2 ranks，OMP／OpenBLAS 為 1，`IGA_PROFILE=1`，KSP GMRES＋LU/MUMPS，
`ksp_rtol=1e-12`。場比較沿用 relative L2 `1e-6`／零參考 absolute L2 `1e-12`。
這些是正確性與失敗處理測試，非效能或跨節點驗收。

第一次 `native/` 與第二次 `native-final/` 選到原本供 adapter 合約使用的
幾何；即使提供非零 profile，該 CLI 路徑仍收斂於零場，沒有 Newton update。
控制器正確拒絕「iteration」注入實際落在 convergence 的結果，兩份 failed
summary 與 logs 保留；沒有把 convergence 的通過當作 iteration 的覆蓋。

`native-driven/` 改用已驗證的 VCA flow geometry，在獨立副本切換 flow_only、
移除 external circuit／coupling、設定非零入口流；另建立 resistance outlet
副本。控制器要求 baseline 確實執行非零次數的線性求解，才進入故障矩陣。
Transport fixture 設定非零 prescribed velocity 與不同初始／入口濃度，確保
時間步有實際演化。原始 fixtures 不修改，副本 hashes 保存在 summary。

健康比較使用本輪封存的 binary，binary hashes 在 summary。重試是新作業，
不宣稱已損壞的同一 stdout writer 可以自行恢復。若失敗發生在最後摘要或
profile，數值輸出可能已經完整寫出；以全部 rank 的非零退出判定 job 失敗。
不撤回已寫入的字元，也不將本輪當作 HPC-05 的原子發布／checkpoint bundle。

## 最終結果

`native-driven/summary.json` 全部通過：412 個 MPI 作業，含 200 個故障、
200 次健康重試、8 個修改前參考及 4 個 restart seed；共核對 618 份 rank
報告。全部故障在指定的共同階段退出 1，全部健康作業退出 0，沒有 timeout。
310 個 velocity／pressure／species 最終場比較的 relative L2 與 absolute L2
皆為 0；數值 gate 保持原值。

Flow 覆蓋 11 個訊息入口、configured transport 7 個、legacy transport 2 個；
每個入口都驗證 1／2 ranks 與五種 buffer failure。實際 target、stage、argv、
每 rank 退出碼／RSS 與場比較均保存在 summary；精簡計數在 `acceptance.json`。
`runtime-results.json` 記錄另行通過的 333 個 trial／port／staged checks。

本批只改變 stdout failure 的偵測與共同退出。每階段立即 flush 會改變輸出
時序，本輪沒有宣稱端到端效能不變；組裝／求解 timing 及 RSS 原始值保留於
rank reports，未將並行正確性測試當作 performance benchmark。

## 重現

```bash
make -C solvers/cpu PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real \
  iga_navier_stokes iga_solve iga_transport \
  flow_stdout_failure_test transport_stdout_failure_test legacy_stdout_failure_test
python3 scripts/hpc_solver_stdout_regression.py \
  --flow-case outputs/hpc01/checkpoint-write/vca/vca_3d_smoke_test/tubularflowiga-vca-3d-smoke \
  --transport-case outputs/hpc01/stream-boundaries/staged-final/1-1 \
  --legacy-case outputs/hpc01/tool-assets/tools/fixture \
  --baseline-dir outputs/hpc01/solver-stdout/before \
  --output-dir /path/to/new-output
make -C solvers/coupling three_d_trial_failure_test three_d_port_failure_test three_d_staged_failure_test
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  PETSC_OPTIONS='-ksp_type gmres -pc_type lu -pc_factor_mat_solver_type mumps -ksp_rtol 1e-12' \
  timeout --kill-after=5s 180s mpiexec --oversubscribe -np 3 \
  solvers/coupling/three_d_trial_failure_test /path/to/new-runtime-output
```

其餘兩個 runtime tests 使用相同環境、對應 binary 與不同新目錄。
Trial／port／staged 的 world 3-rank 及 split 1／2-rank 回歸已全部退出 0，
原有 rollback／retry、port／species 守恆及 serial/group 比較通過。
上述 runtime 合計 333 項既有檢查通過。CPU 與 runtime tests 建置均無 compiler warning。未更改 CUDA 或 immersed／FSI
求解器，未再跑 GPU／完整 FSI；上一批對應驗收見 `HPC_01C_SERIAL_TOOL_PROGRESS.md`。

剩餘：接續 F02 的 1D／graph／sequential／CUDA 入口，以及
尚未覆蓋的串行工具、F05／F06 與 HPC-01D。整份 38 項目標繼續保留。
