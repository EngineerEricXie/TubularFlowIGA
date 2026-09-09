# HPC-01C：graph 啟動的 PETSc 選項一致性

日期：2026-09-08。狀態：**本報告範圍驗收通過；HPC-01C 整體仍部分完成。**
接續 [配置與建構檢查](HPC_01C_CONFIGURATION_PROGRESS.md)。

## 實作與範圍

`CollectivePetscOptions.hpp` 在 supplied communicator 內比較啟動時可見的
PETSc option names、values 與 flags。generic／bifurcation graph 在建立
runtime 前呼叫此檢查；包含 prefixed、未知及尚未使用的選項。
本地擷取錯誤先由 `CollectiveLocalStage` 協調，再進入分塊文字比較。

只排除 `--graph-case`、`--output-dir`、`--stop-after-step`、
`--three-d-max-newton` 與 `-options_file`。前四項由既有 graph 參數／配置
檢查處理；最後一項是已載入的檔案位置，其載入後的選項仍參與比較。
因此相同輸入可位於不同本地目錄，檔案內容造成的 solver 設定差異仍會拒絕。

實作使用公開的 `PetscOptionsGetAll`、`PetscOptionsLeftGet/LeftRestore`、
`PetscOptionsFindPair` API，沒有依賴 PETSc 私有結構。
不能直接比較或重新解析 GetAll 的顯示字串：一個含空白及選項文字的 value，
可能與兩個獨立選項產生完全相同的顯示。helper 以確切 value 長度走訪顯示，
再將 names 與 values 排序、加上長度及 flag 標記；測試明確重現並拒絕此碰撞。
未使用項目由 LeftGet 讀取，FindPair 只查原本已使用的項目，保留 unused tracking。
RAII 釋放 PETSc 提供的 buffers；無法辨識的 display encoding 明確失敗。

此保證是載入後的**選項資料庫快照**，不是所有已建立 KSP/PC 物件、library
state 或稍後修改的等價證明；同值的不同文字表示也不做數值正規化。
後續 runtime 建構的再次檢查已沿用 graph 的應用程式參數排除清單，並以
含 flow／transport 的原生副本 CLI 驗證，見 [資產一致性驗收](HPC_01C_ASSET_PROGRESS.md)。
本機 PETSc 3.15.5 已實測，其他版本仍須驗證顯示格式與 API 行為。
不改物理公式、檔案格式、trial／commit／rollback 或原有數值門檻。

## 命令與環境

WSL 工作站 Intel i9-14900KF，GCC 11.4、Open MPI 4.1.2，
PETSc 3.15.5 real/double、32-bit PetscInt、MUMPS。
HEAD `ee8da2ab528849814b6d12c8186476b5f7f59124` 加未提交工作樹；
來源與執行檔身分以本批 inventory 與 summary hashes 為準。

```bash
make -C solvers/cpu collective_failure_test collective_petsc_options_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
make -C solvers/coupling multidomain_failure_test multidomain_subcommunicator_test petsc \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
python3 scripts/hpc_failure_regression.py --output-dir /tmp/hpc-petsc-options-new-run
```

harness 要求新的輸出目錄，固定 3 ranks、OMP／OpenBLAS threads 各 1，
使用 preonly／LU／MUMPS。child timeout 60 秒、report rendezvous 15 秒、
job timeout 90 秒與 kill grace 5 秒；預期非零退出的原始 rank reports 保持 failed。
subcommunicator 與原生 MPMD CLI 的完整命令另存於下述 JSON。

## 結果

最終 `verified/summary.json` 的 **22 個案例全部通過**：原有 15 項，加上
options-unit 與六個 graph 案例。

| 新案例 | 結果 |
|---|---|
| options-unit | world 與 1+2 groups 驗證插入順序、used flags、空值、特殊字串、顯示碰撞及跨 rank 差異 |
| petsc-ksp、petsc-prefix、petsc-flag | rank 1 設定不同，三 rank 均回傳 1，MPI 作業退出 1 |
| petsc-used | 相同選項、不同 used metadata，正常求解且保留 unused tracking |
| petsc-file | 三個不同位置的相同 options files，正常求解 |
| petsc-file-different | rank 1 檔案載入不同 rtol，三 rank 均回傳 1 |

故障案例皆在建構 solver 前以共同診斷退出，沒有 timeout 或成功 manifest。
整批最大單 rank wall `0.714781 s`、peak RSS `45,555,712 bytes`；
這是微型正確性測試，包含啟動成本，並非 scaling 或獨立效能重複量測。

另以真正的 `iga_multidomain_flow` MPMD 啟動驗證 PetscInitialize 的 argv 路徑：
不同 graph 絕對路徑的相同副本可完成，四個 edge pressure／flow 欄位相對誤差
均為 0；僅 rank 1 使用 `-ksp_type cg` 時作業退出 1，無成功 manifest。

`multidomain_subcommunicator_test` 的 1+2 群組刻意使用不同 `-ksp_rtol`。
flow 與 species 皆通過，兩 rank 群組最大相對 L2 分別為 `7.45082e-15`、
`3.35727e-15`，單 rank 群組為 0。既有 CPU `1e-6` 門檻及物理 gate 不變。
受影響 binary 重建無新增 compiler warning。未重跑 CUDA、1005-node body case
或完整 multidomain CLI suite；前一報告的這些結果不作為本版本的新執行證據。

## 證據與後續

本地 ignored 產物位於 `outputs/hpc01/petsc-options/`：
`verified/summary.json`、每 rank 的 reports／logs、`subcommunicator.json`、
`native-cli.json`、建置 logs、`evidence-summary.json` 與 `inventory.json`。
`initial` 保留最終 unit test 補強前的結果；權威案例集合為 `verified`。
已重算並核對 summary 所記 binary hashes、回歸 log hashes 及預期退出碼。

HPC-01C 剩餘其他 CLI、外部資產內容身分、後段建構與耦合／輸出失敗邊界，
仍依 [完整待辦清單](../WORKSTATION_HPC_TODO.md) 推進。本協議要求全員同序
呼叫，不能將含 collective 的大函式放入 local callback，也不承諾 process-loss 恢復。

API 參考：[GetAll](https://petsc.org/release/manualpages/Sys/PetscOptionsGetAll/)、
[LeftGet](https://petsc.org/release/manualpages/Sys/PetscOptionsLeftGet/)、
[LeftRestore](https://petsc.org/release/manualpages/Sys/PetscOptionsLeftRestore/)、
[FindPair](https://petsc.org/release/manualpages/Sys/PetscOptionsFindPair/)。
上述線上 release 文件與本機 3.15.5 不同；版本相容性以本機建置及執行證據為準。
