# HPC-01C：1D checkpoint／restart 一致失敗處理

日期：2026-09-08。狀態：**本批 25 項 CLI／格式案例及單／三 rank unit 通過；
HPC-01C 與 HPC-05 整體仍未完成。**

## 實作與相容性

`OneDCheckpoint.hpp` 現在先協調本地 metadata 讀取、解析、layout 與 index
容量檢查，再比較同一 communicator 的 metadata 快照。狀態讀入後，檢查大小、
有限值與正面積，確認各 rank 的完整狀態副本一致，最後協調 unpack。
計算 packed size 時使用有界的 unsigned 運算，避免原本 `3*cells` 等 signed
integer overflow。CLI 的 restart 狀態複製、驗證／restore，以及 checkpoint
目錄與 metadata 建立也已加入本地錯誤協調。

既有 1D runtime 已在每 rank 保留完整狀態。讀取時各 rank 以 `COMM_SELF`
載入自己的 checkpoint 副本，寫出時由 communicator rank 0 以 `COMM_SELF`
寫入完整向量。PETSc 本地 I/O 呼叫結束、清理完成後，才透過 supplied
communicator 協調結果；沒有把 group collective 藏進 local callback。
這讓缺檔、截斷或 viewer 錯誤不會把其他 rank 留在全域 viewer collective 裡。

PETSc 呼叫逐一檢查回傳碼，局部使用 `PetscReturnErrorHandler` 並在離開時
還原 handler。RAII 管理部分建立的 Vec／viewer 及 array read view；正常路徑
明確檢查 close，metadata 也在 close 後確認寫入成功。
兩個 public checkpoint 函式新增可選的 borrowed communicator，舊參數形式
保留 `PETSC_COMM_WORLD` 預設；函式不取得或釋放 caller 的 communicator。

metadata schema 2、packed 順序及 PETSc binary `.state` 格式保留。
checkpoint viewer 固定 filename／header 契約，略過 viewer options 及 `.info`
內的選項，避免由附屬檔案改變載入行為。沿用以 PREFIX 定位 `.json`／`.state`
的既有介面；此批沒有重新定義 metadata 的 `state_file` 欄位。

此策略仍有完整狀態複本與重複讀取成本，寫出時 root 另持有完整 sequential Vec。
本批小型案例不證明這適合大型 1D network；HPC-06 仍需量測及改善其記憶體／I/O。
相同 schema 的全員一致但有限值損毀，沒有 checksum 可供辨識；狀態與 metadata
也尚非原子發布，不保證寫入失敗後可保留上一份 checkpoint。這些是 HPC-05
版本化／發布協議的待辦，不因本次錯誤退出驗證而視為完成。

## 命令與環境

```bash
make -C solvers/one_d iga_1d one_d_petsc_test one_d_checkpoint_format_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
python3 scripts/hpc_one_d_checkpoint_regression.py --output-dir /tmp/hpc-checkpoint-new-run
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 timeout 60 \
  mpiexec --oversubscribe -np 1 solvers/one_d/one_d_petsc_test
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 timeout 60 \
  mpiexec --oversubscribe -np 3 solvers/one_d/one_d_petsc_test
```

unit 的兩次執行要依序進行，其既有 fixture 使用固定暫存檔名。
WSL、Intel i9-14900KF、GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5
real/double、32-bit PetscInt、MUMPS；CLI harness 使用 3 ranks，OMP／OpenBLAS
threads 各 1，preonly／LU／MUMPS。HEAD
`ee8da2ab528849814b6d12c8186476b5f7f59124` 加未提交工作樹，來源身分見 inventory。
raw binary 故障注入會先檢查 real/double、32-bit 格式，不將此測試宣稱為其他
PETSc scalar/index 配置的驗收。重建完成，無新增 compiler warning。

## 結果

`clock-verified/summary.json` 的 **25 項全部通過**：

| 案例 | 驗收 |
|---|---|
| rigid full／save／resume | 原生 10 步，中途第 5 步存檔，續跑到第 10 步 |
| explicit full／save／resume | 原生 compliant bifurcation 20 步，第 10 步存檔續跑 |
| species full／save／resume | 原生六物種、PETSc pressure network 120 步，第 60 步存檔續跑 |
| legacy-format、legacy-resume | 新檔經原本的 distributed VecLoad／VecView 讀寫後逐位元相同，再載入續跑 |
| metadata missing／malformed／different／overflow | rank 1 缺檔、JSON 截斷、合法但不同 inlet flow、超過 index capacity |
| state missing／truncated／header／nonfinite／area／different | rank 1 缺檔、資料截斷、錯誤 class ID、NaN、零面積、合法但不同面積 |
| write parent／state／metadata／buffered | rank 1 目錄不可建立；root state／metadata 位置是目錄；metadata 指向 `/dev/full` |

正常與續跑的 pressure、flow、area、species、derived 欄位及幾何／場身分一致，
**所有非時間欄位誤差為 0**。剛體的時間也完全相同；explicit／species 的時間
欄位最大相對誤差分別 `4.60348e-17`／`6.78149e-17`。
regular output 使用 step*dt，restart 初始輸出使用累積 physical time；比較沿用
CLI `ValidateRestart` 的 `1e-12*max(1, |t1|, |t2|)` clock gate。
CPU 場 relative L2 `1e-6`、零參考 absolute L2 `1e-12` 均未放寬。

初輪 `initial` 因 exact timestamp 篩選漏掉 checkpoint 那一列而停止，
不是場數值 gate 通過；此紀錄保留。對齊既有 restart clock 規則後完整重跑，
以 `clock-verified` 為權威結果，不抹去初輪失敗。

所有 14 項故障的三份 rank reports 與 MPI launcher 均退出 1，沒有 timeout，
root 診斷符合指定階段，沒有成功的 simulation summary。正常案例退出 0，
驗證停止／恢復步數及完整剩餘 trajectory；比較物種和 derived CSV 的名稱與順序。
child timeout 60 秒、job timeout 90 秒、kill grace 5 秒；wrapper 保留既有
15 秒 report rendezvous 與原始非零報告。
最大單 rank wall `0.414247 s`、peak RSS `40,951,808 bytes`，僅為小案例
端到端量測，不是效能重複統計或 scaling 證據。

既有 `one_d_petsc_test` 單／三 rank 均通過，保留四種 implicit formulation
與 pack/unpack 檢查。新增的 1+2 communicator 群組以不同壓力、不同呼叫次數
獨立寫入／讀回 checkpoint，packed 值逐位元相同；兩 rank 群組也確認 local
rank 1 缺 metadata 時，組內全員收到錯誤，不影響另一群組。

## 證據與剩餘範圍

本地 ignored `outputs/hpc01/one-d-checkpoint/` 保存 `initial`、
`clock-verified/summary.json`、每 rank reports／logs、輸出場、checkpoint、
legacy `.state`、`unit-groups-summary.json`、build logs、`evidence-summary.json`
及 `inventory.json`。已核對目前 production／format／unit binary hashes、
預期退出碼及 log hashes；先前 unit（未加入群組測試）的紀錄另行保留。

後續已補上 CLI 時間步與 combined trial 的局部協調，見
[trial 進度](HPC_01C_ONE_D_TRIAL_PROGRESS.md)。HPC-01C 仍需 implicit runtime
內部、staged 1D 與其他 CPU／coupling CLI 的失敗邊界；HPC-05 仍需跨 domain 狀態契約、checksum、
原子發布、中斷恢復及新 runtime 整合。本批沒有 CUDA 或跨節點驗收。

API 依據：[PetscViewerBinaryOpen 的 collective 契約](https://petsc.org/release/manualpages/Viewer/PetscViewerBinaryOpen/)、
[PetscReturnErrorHandler](https://petsc.org/release/manualpages/Sys/PetscReturnErrorHandler/)。
線上 release 文件版本較新，實際相容性仍以上述本機 3.15.5 的建置與測試為準。
