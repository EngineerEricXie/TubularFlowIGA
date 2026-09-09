# HPC-01C：文字輸出的關檔錯誤

日期：2026-09-08。狀態：本批驗收通過；HPC-01C 與整份清單仍未完成。

## 問題與修正

`WriteVtu`、`WritePvd`、`WriteVelocityManifest` 原先只在 stream 尚未關閉時
檢查狀態；小檔案可能仍在緩衝區，關檔才發生的寫入失敗會被漏報。
`WritePhysiologyManifest` 則連寫入後的狀態都未檢查。

四個共用 writer 現在明確 `close()`，再檢查 stream 狀態並拋出包含輸出路徑的
錯誤。資料、欄位順序、數值精度與檔案格式不變；呼叫端既有錯誤協調可接收到
這些例外。沒有更動方程、求解器或數值驗收門檻。

這不等於原子發布或斷電持久化：失敗後可能保留截斷檔案，且沒有新增 `fsync`。
本批亦未完成所有呼叫端的 MPI 錯誤協調，不能推論每一個 CLI 都已可靠退出。

## 故障重現與驗收

新增 `solvers/cpu/tests/test_text_output_failure.cpp`，接入 CPU `make test`。
每個 writer 分別測試：小檔案在關檔時失敗、大檔案在串流寫入時失敗，以及
輸出目標是目錄而無法開啟。前兩者以子程序設定 `RLIMIT_FSIZE=64` 並忽略
`SIGXFSZ`，實際產生 64-byte 的截斷 regular file，要求 writer 回報指定錯誤。
限制只套用子程序，測試不需填滿磁碟，也不改使用者的程序限制。

- 修改前的四個 writer 共漏報 5 個案例：4 個關檔錯誤與 physiology 串流錯誤。
  保存的舊標頭編譯版本實際退出 1，並保留失敗輸出。
- 修改後 12 個故障案例與 8 次重新寫入全部通過。重試使用原先截斷檔案，
  完整內容與健康輸出逐位元組相同。
- 四種健康輸出與修改前逐位元組相同，並另用 XML／JSON／CSV reader 解析。
- 既有 `vtk_output_test` 通過。新測試使用明確條件檢查，不依賴可被 `NDEBUG`
  關閉的 assert；單次測試有 30 秒 timeout。非 POSIX 平台回報 SKIP 並退出 77，
  不冒充完成故障注入驗收。

## 呼叫端回歸

| 驗收 | 結果與範圍 |
|---|---|
| CPU configured transport | 81 筆觀察、154 份 rank reports，67 份場比較通過；其中 12 份是參考自比 |
| CPU flow | 51 筆觀察、94 份 rank reports，32 份指定場比較通過；保留 legacy 跨 rank 壓力的 4 筆未通過觀察 |
| CUDA transport | 本機 RTX 4080 SUPER，prescribed／非恆定 snapshot series／VTKHDF 三個模式，各比較初始、兩步及最終文字場，共 12 份 CPU/GPU 場比較通過 |

CPU transport 最大 relative L2 為 `3.0631271369706667e-12`；flow 的指定比較
最大值為 `1.8156697926380728e-13`。CPU 門檻仍為 relative `1e-6`，零參考
absolute `1e-12`。legacy 的近零壓力跨 rank 問題仍保留，未以同 rank 相容性
取代其失敗證據，見 [flow 輸入報告](HPC_01C_FLOW_INPUT_PROGRESS.md)。

CUDA 最大 relative L2 為 `2.7219049804952843e-6`，通過既有 `1e-5` 門檻；
field-name sidecar 一致。這裡的 VTKHDF 模式比較指文字場，不是新做的 HDF5
dataset 比較或 ParaView 驗收。本批沒有修改 HDF5 writer。

CPU 的輸入、binary 與原始 rank 日誌雜湊均重新核對；健康／預期失敗退出碼與
timeout 狀態符合各案例條件。這些小案例部分重疊執行，只用於正確性回歸，
不提供效能、OpenMP 加速或跨節點擴展性證據。

## 重現命令與證據

本機環境：GCC 11.4、OpenMPI 4.1.2、PETSc 3.15.5 real64／Int32、MUMPS；
CUDA 12.6、SM 89。CPU CLI 的 OMP／BLAS 固定 1，使用 core mapping／binding。
CPU flow 的 OpenMP 執行檔亦重建，但本批未新增該版本的執行緒驗收。

```bash
make -C solvers/cpu text-output-failure-test vtk_output_test
./solvers/cpu/vtk_output_test
make -C solvers/cpu iga_solve iga_navier_stokes iga_navier_stokes_openmp \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
python3 scripts/hpc_transport_cli_regression.py \
  --case-dir outputs/hpc01/checkpoint-write/verified/staged-failure/1-1 \
  --reference-binary outputs/hpc01/text-output-close/iga_solve-before \
  --output-dir outputs/hpc01/text-output-close/transport
python3 scripts/hpc_flow_input_regression.py \
  --case-dir outputs/hpc01/assets/vca-cli/tubularflowiga-vca-3d-smoke \
  --reference-binary outputs/hpc01/text-output-close/iga_navier_stokes-before \
  --output-dir outputs/hpc01/text-output-close/flow
conda run -n tubularflow-cuda make cuda CUDA_ARCHS=89
```

回歸目錄已存在，重跑須指定新目錄。完整 CUDA 命令、退出碼與場比較位於
`outputs/hpc01/text-output-close/cuda-result.json`；其餘證據包括
`before-four-result.json`、`after-result.json`、`audit.json`、兩組 CLI 的
`summary.json`、build logs 與 `acceptance.json`。正式來源與執行檔另行封存，
不提交測試產物。

## 剩餘工作

繼續補上 flow 建構後的本地配置／時間步／輸出邊界，以及其他 CLI、adapter
和 executor 的剩餘協調。共用 HDF5、geometry report、其他文字輸出的完成
檢查亦須逐一核對。整體持久化協議與分散式輸出仍分別由 HPC-05／06 驗收。
