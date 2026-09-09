# Native graph checkpoint／restart

`iga_multidomain_flow` 與 `iga_1d_3d_bifurcation` 已接入 accepted-step checkpoint。
目前支援原生 0D／1D／貼體 3D flow graph，以及 1D／貼體 3D species graph；
恢復使用相同 rank 數與 communicator membership。浸入式、移動、FSI、CUDA graph
及重分區恢復由 HPC-05D 追蹤。單機驗收見
[native graph 報告](progress/HPC_05C_NATIVE_GRAPH_PROGRESS.md)。

## 保存與新的作業續跑

使用同一份輸入、相同執行緒及 PETSc 數值選項，且各 `.ntiga` 分區符合 rank 數。
Makefile 在編譯時嵌入來源內容 SHA，包含尚未 commit 的修改，不依賴執行時原始碼。
自訂建置須同樣設定 `IGA_NATIVE_CHECKPOINT_SOURCE_SHA256`，可用
`python3 scripts/native_checkpoint_source_identity.py` 取得；缺少指紋的程式仍可正常
執行未啟用 checkpoint 的工作，但會拒絕保存與載入。

```bash
make -C solvers/coupling iga_multidomain_flow PETSC_DIR="$PETSC_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps'

# 首次作業：每三個 accepted macro-steps 保存，並於第六步停止。
mpiexec -np 3 solvers/coupling/iga_multidomain_flow \
  --graph-case CASE --output-dir FIRST_OUTPUT \
  --checkpoint-dir CHECKPOINT_ROOT --checkpoint-every 3 --stop-after-step 6

# 新作業：恢復最後完整世代，完成原設定剩餘步數。
mpiexec -np 3 solvers/coupling/iga_multidomain_flow \
  --graph-case CASE --output-dir RESUMED_OUTPUT \
  --restart-dir CHECKPOINT_ROOT --checkpoint-dir CHECKPOINT_ROOT --checkpoint-every 3
```

`FIRST_OUTPUT`／`RESUMED_OUTPUT` 必須尚不存在。總時間步設定仍來自原 case；
`--stop-after-step` 只縮短本次作業的 horizon，不能短於恢復世代。
新的輸出包含 checkpoint 保存的完整 accepted history prefix，再接新步數；不依賴
先前輸出目錄，也不重複已接受的步數。若只提供 `--restart-dir`，來源 bundle 保持唯讀。

| 參數 | 行為 |
|---|---|
| `--checkpoint-dir ROOT` | 保存至共享的 POSIX 目錄；可建立缺少的目錄 |
| `--checkpoint-every N` | 正整數，預設 1；必須同時提供保存目錄 |
| `--restart-dir ROOT` | 載入最後完整且相容的世代 |
| `--stop-after-step N` | 到指定 accepted step 停止；該最後步即使未落在頻率上也保存 |

正常到達 case 最後步也會保存。每個世代保存完整 graph state、下一步 pressure
guess、species donor hysteresis 與 accepted history；owned 3D fields 由各 rank 寫入，
小型 replicated state／history 由 group root 寫入。Aitken 的 step-local history 依原演算法
重新初始化。時鐘使用原始逐步累加值。

## 排程預警

提供 `--checkpoint-dir` 時，程式從輸入載入前開始接收 `SIGUSR1`。handler 只設定旗標；
下一個成功的 accepted macro-step 完成後，全群保存並正常退出，即使未到保存頻率。
任何一個參與 rank 收到訊號都會觸發這個行為。預警須送到實際 solver rank，並留足
一個完整 macro-step 加保存時間；若下一步失敗或遭強制終止，只能恢復先前完整世代。
launcher／batch shell 的訊號轉送依 scheduler 配置而異，不能假定只通知 shell 即有效。
實際站點整合與多節點排程仍依 [HPC-09](WORKSTATION_HPC_TODO.md) 驗收。

## 相容性與失敗處理

輸入內容 SHA 包含 graph／domain configuration 和所有外部資產；case 可整體搬移，
內容必須相同。另核對來源 SHA、MPI／PETSc 版本、scalar／index 寬度、communicator
成員與 rank 數、compiler／OpenMP／fast-math 身分、執行緒環境、PETSc 數值選項與
native Newton／守恆門檻。來源指紋保守納入共用 header 集合；即使修改不影響此 case，
仍需重新驗證相容性。保留同一建置與數值配置，第一版不承諾跨 build 自動轉換。

先寫分片、校驗與 `fsync`，最後發布 manifest。失敗世代留在原位；載入跳過不完整、
截斷、損毀或缺片世代，回到最後完整且相容的 checkpoint。同一步重新保存會選擇
新的 epoch ID，不覆蓋既有檔案。相容性或 candidate 驗證失敗時，不執行新步、不發布
graph 輸出。每份 graph 候選失敗後須丟棄，這不是已運行 graph 的原子替換介面。

格式與耐久性條件見 [bundle v1](architecture/COUPLED_CHECKPOINT_BUNDLE.md)。目前為
共享 POSIX 檔案系統；尚未宣稱大型／跨節點檔案系統的數值、耐久性或吞吐驗收。
history 按 record 串流，每筆 metadata 上限 16 MiB，完整 prefix 可超過該上限；
既有記憶體內 history 與每代完整 prefix 的成長成本仍由 HPC-06 處理。
