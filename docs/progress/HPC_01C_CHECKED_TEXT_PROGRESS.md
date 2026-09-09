# HPC-01C：完整文字讀取與序列化錯誤

日期：2026-09-08（EDT）。本批指定範圍通過；其餘共用串流與入口稽核仍未完成。

## 問題與修正

`ostringstream << input.rdbuf()` 的複製錯誤會反映在目的串流上，輸入的
`good()` 仍可能為 true。若檔案在完整 JSON 加尾端空白後發生讀取錯誤，
舊式檢查可能把已讀取的有效前綴交給 parser，而未拒絕不完整的讀取。

新增 [CheckedText.hpp](../../include/CheckedText.hpp)，以分塊 `istream::read`
讀取並檢查 bad／fail／EOF；字串 append 的配置失敗直接拋出。保留空輸入、
最後不足一塊的正常 EOF、嵌入 NUL 的原始位元組與各 caller 的 open/path
政策。Helper 本身不做 MPI，既有 caller 的本地共同階段負責協調錯誤。

已替換 17 個讀取點：case／simulation／mesh／1D／multidomain／0D 配置、
浸入式 JSON、temporal CSV、velocity manifest、coupling replay、flow／
transport／VCA checkpoint metadata、Womersley reference，以及原生 1D、
sequential 與 graph 的 ReadText。Config checker 的檔案與 stdin 入口也使用
同一 helper；沒有更改 parser、schema 或缺檔時的既有 optional 行為。

以下中間串流啟用 badbit／failbit exceptions，避免默默回傳截斷的字串：

- CPU flow／transport 執行控制值。
- Sequential／multidomain JSON escape，sequential strong CSV header／row。
- Flow／transport／VCA checkpoint metadata、VCA device identity／fingerprint
  與 time-indexed output suffix。

VCA fingerprint 的逐字元讀取亦檢查正常 EOF。既有 hash 演算法、浮點格式、
CSV／JSON schema 與檔案命名保持不變。

## 驗收

證據根目錄為 `outputs/hpc01/checked-text/`。

| 範圍 | 實際結果 |
|---|---|
| Checked reader | 7 個 byte／EOF 邊界案例；3 種舊式複製漏報重現；標準例外、bad_alloc、非標準例外在讀取前／讀取有效前綴後共 6 次拒絕與健康重試；預先 fail 的 stream 亦拒絕 |
| 原生配置失敗 | 直接呼叫 sequential 的 ReadText／JsonEscape 與三種 checkpoint serializer，僅最後一個 rank 注入首次配置失敗；world 3 ranks 與獨立 1／2-rank 子群，共 15 次共同拒絕與 15 次重試 |
| CPU／1D／coupling／mesh | CPU `test`、1D `core-test`、coupling `test`、mesh／mesh-test 通過；後續 checkpoint serializer 修改後，再通過三個 metadata unit tests |
| Sequential | 完整 MPI smoke 通過，含 explicit／graph、fixed／Aitken 與 subcycling |
| Configured transport CLI | 81 個案例、154 份 rank reports 通過，涵蓋讀取、輸出、續跑、錯誤拒絕及既有場比較 |
| Flow CLI | 51 個案例、94 份 rank reports 符合既有判定，包含原有 legacy 跨 rank 壓力差異的保留觀察；不宣稱那些既有數值限制已修復 |
| Multidomain 原生相容性 | flow、species、0D graph 各跑 1／2／3 ranks，新舊共 18 次執行，57 份 CSV 逐位元組相同 |
| CUDA | CUDA 12.6／sm_89 建置通過；RTX 4080 SUPER 的新舊各四個 transport 場相同，最終 CPU/GPU relative L2 `2.397514096418484e-6`，通過既有 CUDA `1e-5` 門檻 |

配置故障測試的注入只在本地文字操作內開啟；MPI／PETSc 初始化與錯誤協調
不在注入區間。讀取故障使用會拋出例外的 streambuf；這不是所有 OS read
呼叫點或所有配置點的窮盡故障測試。

首輪 CUDA 比較腳本誤用了 CPU 跨 rank 的 `1e-6` 門檻，返回 1；新舊 GPU
的誤差與場檔完全相同。既有
[基準 catalog](../../benchmarks/hpc_baselines.json)與
[HPC-00D 報告](HPC_00D_REPORT.md)原先即規定 CPU/CUDA 為 `1e-5`。
保留 `cuda-run.log`、原始輸出與較嚴格比較，最終重新建置後的結果在
`cuda-final-runs.json`／`cuda-final-comparison.json`。未改 catalog、solver
停止條件、物理模型或原有驗收門檻。

另保留兩個測試準備失敗：flow controller 首次使用不符合命名規則的 fixture，
在啟動 solver 前因缺少 `fixture.ntiga` 退出；擴充 native test 時漏 include
兩個 checkpoint headers，已補齊後重建。最終 production／unit builds
無新增 compiler warning。

## 重現命令

本機 WSL、GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real64／Int32。
CPU native comparison 使用 GMRES＋MUMPS LU、`ksp_rtol=1e-12`；各既有
controller 沿用自己的已記錄設定，OMP／BLAS 均固定 1。這些是相容性驗收，
並非隔離效能量測；沒有新增跨節點或 FSI 分散計算宣告。

```bash
make -C solvers/cpu -j2 test iga_solve iga_navier_stokes \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
make -C solvers/one_d -j2 iga_1d core-test
make -C solvers/coupling -j2 iga_1d_3d_explicit iga_multidomain_flow \
  iga_1d_3d_bifurcation sequential_text_failure_test test
make mesh mesh-test
conda run -n tubularflow-cuda make cuda CUDA_ARCHS=89
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  timeout --kill-after=5s 60s mpiexec --map-by core --bind-to core -np 3 \
  solvers/coupling/sequential_text_failure_test /absolute/path/to/new-output
```

原生 CLI 矩陣的 controller、完整 argv、輸入 hashes 與每 rank 狀態保存於
`transport-cli/summary.json`、`flow-cli-final/summary.json`；graph 命令與比較
在 `graph-runs.json`／`graph-comparisons.json`。重跑需指定新的 output directories。
CUDA 的原始與修正後比較腳本亦保存在證據根目錄。

## 後續追蹤

- 本批結束時尚未補強的 `OneDCheckpoint.hpp` network fingerprint、metadata
  serialization 與 `rdbuf()` metadata 讀取，以及 VTK signature／檔名與資源摘要，
  已在後續 [文字 helper 驗收](HPC_01C_TEXT_HELPER_PROGRESS.md) 完成指定修正與測試。
  這些後續結果不屬於本報告的原始 source archive。
- 逐一核對剩餘 helper 的 caller 是否在適當的本地共同階段；本地 callback
  不得包住 MPI collective，也不能以 grep 沒有匹配代替完整入口稽核。

HPC-01C／01D 及後續階段仍依原清單，整體保持 9／38 項勾選。
