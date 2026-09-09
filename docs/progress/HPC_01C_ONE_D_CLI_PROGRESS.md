# HPC-01C：原生 1D CLI 的輸入與輸出失敗邊界

日期：2026-09-08。狀態：**本報告列出的 18 項案例通過；HPC-01C 仍部分完成。**

## 實作

`iga_1d` 在進入求解前，協調參數解析、配置讀取、network／runtime 建構及
initial state 的本地例外。各 rank 的 system 選擇、check 模式、停止步數、
checkpoint 頻率與 restart／checkpoint 是否啟用必須一致，避免部分 rank
提早 finalize 或走進不同的 collective 分支。

配置比較使用實際供 parser 解析的 JSON 快照；格式差異也視為不同輸入。
相同副本可放在不同本地目錄。啟動時也比較可見 PETSc option entries，
沿用 [graph 的選項擷取協議](HPC_01C_PETSC_OPTIONS_PROGRESS.md)，
排除由應用程式處理的七個 `--` 選項及已載入的 `-options_file` 位置。
其餘選項包含 flags；沒有聲稱比較 PETSc 物件的所有有效預設或後續修改。

root 的 writer 建立、初始場、時間步場與最後結果寫出均由所有 rank
共同確認後才前進。`OneDOutputWriter` 現在在每個輸出步 flush CSV 並檢查
結果；VTP、PVD 與 summary 在 close 後檢查狀態，CSV 的 close 也在 summary
寫出之前確認。這能偵測原本只在 buffer flush 時發生、沒有被回報的寫入錯誤。
未加入 fsync、原子發布或完整的持久化協議；這些仍由 checkpoint／I/O 階段追蹤。

local callbacks 只含本地工作。runtime 建構只是儲存 implicit advance callback，
initialization 使用既有本地 rigid 初始化；`SolveTrial`、PETSc checkpoint
read/write 仍在 callback 外，沒有將含 collective 的大函式包成 local stage。
物理公式、schema、檔案內容格式與數值門檻不變。flush 與協調會增加 I/O／通訊
成本，未宣稱速度改善；後續效能矩陣必須包含此成本。

## 驗收命令與環境

```bash
make -C solvers/one_d iga_1d \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
python3 scripts/hpc_one_d_cli_regression.py --output-dir /tmp/hpc-one-d-new-run
```

輸出目錄必須不存在。WSL、Intel i9-14900KF、GCC 11.4、Open MPI 4.1.2、
PETSc 3.15.5 real/double、32-bit PetscInt、MUMPS。每 rank OMP／OpenBLAS
threads 各 1，PETSc options 為 preonly／LU／MUMPS。
HEAD `ee8da2ab528849814b6d12c8186476b5f7f59124` 加未提交工作樹，精確來源與
binary hashes 見 inventory。編譯通過，沒有新增 compiler warning。

測試以 MPMD 啟動真正的 CLI，每 rank 使用不同輸入目錄，不在 production
加入測試 knob。child timeout 60 秒、job timeout 90 秒、kill grace 5 秒；
三 rank wrapper 沿用 15 秒 report rendezvous，保留每個 child 的真實非零狀態。

| 案例 | 驗收 |
|---|---|
| rigid-1／3 | 原生 10-step straight tube，單／三 rank 場一致 |
| explicit-1／3 | 原生 20-step compliant bifurcation，單／三 rank 場一致 |
| implicit-1／3 | bifurcation 改用 PETSc pressure_network，dt 0.001、2 steps，單／三 rank 場一致 |
| arguments | rank 1 的非法停止步數，全員退出 1 |
| check-control、stop-control | rank 1 的合法但不同執行分支，全員退出 1 |
| config-missing、config-different | rank 1 缺配置或讀到不同 dt，全員退出 1 |
| geometry-missing | rank 1 缺 SWC，runtime input 階段全員退出 1 |
| petsc-flag | 只有 rank 1 提供 `-ksp_monitor`，選項比較全員退出 1 |
| output-setup | root 的輸出目錄位置是檔案，全員退出 1 |
| output-initial、output-step | 第 0／1 步 VTP 位置是目錄，全員退出 1 |
| output-final | PVD 位置是目錄，在寫出 summary 前全員退出 1 |
| output-buffered | flow CSV 指向 `/dev/full`，初始輸出的 flush 失敗，全員退出 1 |

共 **18 項通過**，各 fault case 的三份 rank reports 均為退出 1、無 timeout，
MPI launcher 也退出 1。檢查指定 root 階段診斷、log hashes 與無成功 summary。
正常案例檢查步數、converged 標記、CSV 列數、有限值、正面積與場身分。
`profile_1d.csv`／`branch_timeseries.csv` 各欄位的單／三 rank 誤差均為 0；
既有 CPU relative L2 `1e-6`、零參考 absolute L2 `1e-12` 不變。
這是平行數值回歸，不是對三種物理模型的新增獨立解析解驗證。

最大單 rank wall `0.314404 s`、peak RSS `40,497,152 bytes`。
這些小案例包含啟動與輸出成本，既不是獨立效能重複測試，也不是分散式加速證據。
rigid／explicit 物理更新仍複製執行，不能將三 rank CLI 可用稱為三倍加速。

## 證據與剩餘工作

產物在 ignored `outputs/hpc01/one-d-cli/`：`verified/summary.json` 包含
完整 argv、輸入／執行檔 hashes、逐 rank 命令與量測；另保留 build logs 和
`inventory.json`。原始 `initial` 在 PETSc 差異案例停止：只改 rank 1 的
`PETSC_OPTIONS` 環境變數沒有得到預期的非零退出，不算該 gate 通過。
最後改用逐 rank 命令列 flag，完整重跑 18 項，以 `verified` 為權威結果。

後續已補上 checkpoint metadata／state 一致性、本地 PETSc I/O 與 restore
錯誤邊界，並驗證完整 multispecies 續跑，見
[checkpoint 報告](HPC_01C_ONE_D_CHECKPOINT_PROGRESS.md)。本報告保留當時 18 項結果。

仍待補齊：合法但不同的 SWC／waveform／replay 內容身分，時間步 inlet／bookkeeping、
`SolveTrial` 與 implicit runtime 內部的階段協調，以及 CPU／explicit coupling CLI。
成功寫出部分檔案也不等於整個目錄原子發布。未在本批重跑完整 1D multispecies、
closed-loop 或 restart suite；這些路徑的完整驗收仍需另行執行。
完整目標依 [待辦清單](../WORKSTATION_HPC_TODO.md) 持續推進。

後續的 CLI 時間步與 combined trial 局部階段已補強並重跑數值回歸，見
[trial 進度](HPC_01C_ONE_D_TRIAL_PROGRESS.md)。implicit PETSc 內部及 staged
1D 邊界仍不可視為完成。
