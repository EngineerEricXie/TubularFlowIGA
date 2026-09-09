# HPC-02 貼體流場 MPI／OpenMP 組裝進度

HPC-02A/B/C 保持未勾選，等待整個階段的 FSI 單 thread 與正式量測驗收。
貼體流場的執行配置、故障回復、獨立 communicator 回歸，以及固定
總核心數的正式重複量測已完成；以下記錄此範圍的證據。

## 實作

`TransientFlowRuntime` 使用既有 `ParallelElementBatch` 與
`ElementAssemblyExecution`，主執行緒先複製 required 節點的本步及歷史場，
worker 計算私有的 Navier–Stokes 元素系統，再由主執行緒按原順序插入
owned rows。保留既有 required 元素覆蓋；不能直接按 METIS owner 篩除，
否則可能遺漏 owned rows 所需的相鄰元素貢獻。壁面牽引仍由主執行緒積分。

執行資源物件在共同錯誤邊界內建構，避免一個 rank 的錯誤設定讓其他
rank 留在 collective。worker 例外在批次全部 join 後，透過既有共同錯誤
協議傳播。失敗後 retry 以 `MAT_FLUSH_ASSEMBLY` 清除未完成的插入模式，
再歸零重新組裝；第一次組裝未填入的預配置位置得以保留。

新增 `iga_navier_stokes_openmp` 可選 CLI；普通 CLI 繼續可以不含 OpenMP。
新增 `parallel_body_fitted_flow_test` 為原生 MPI 故障／場比較測試。
`.ntiga`、外部節點 ID、物理參數及 trial／commit／rollback 介面不變。

## 原生驗證

`outputs/hpc02/hybrid/native-2rank-restored/` 使用直管原生 1,005-node
案例，兩 ranks、各 1／2 threads，LU/MUMPS 與既有 `1e-6` 相對場門檻。
兩 rank 均退出 0；非線性、continuity 與質量 gate 通過。

- rank 1 錯誤的 thread 設定共同拒絕。
- 第一次矩陣組裝於第二批元素 10／12 注入錯誤，批次全部完成後選擇元素 10。
- 未解出的 trial 不可 prepare；rollback 保持已提交狀態，retry 可成功求解。
- 完整稀疏模式建立後再注入錯誤，Abort／新一步重試保持正確狀態。
- 真實 worker 參與、有界暫存及 required 元素覆蓋經過檢查。
- 1／2 threads 的速度相對 L2 為 `2.07708e-19`，壓力為 `2.48147e-18`。

該輪使用封存後恢復的候選來源；已與 root 整合，`native-4rank-root/`
的四 ranks、各 1／2 threads 亦全部退出 0，速度及壓力比較均為 0。
原生斷言與物理收斂 gate 均通過。`native-evidence.json` 與
`root-source-and-binaries.tar.gz` 保留數值紀錄及整合來源／執行檔。
較早的
`native-2rank-candidate/` 兩 rank 退出 195、沒有完成診斷；當時工作回合
曾中斷，僅能記錄未完成，不能以此判定為數值成功或 MPI deadlock。

`subcommunicator-root/` 的既有非 OpenMP 3-rank 回歸也全部退出 0。
1／2-rank 獨立 communicator 各自與單 rank 參考比較，涵蓋流場、傳輸、
積分質量、來源、halo 更新與 checkpoint／rollback；group 0 各場為 0，
group 1 的速度 `3.48315e-15`、壓力 `9.94533e-16`、傳輸 `1.70503e-14`、
質量 `1.62206e-14`，皆通過原有 `1e-6` 門檻。

## 固定總核心數量測

`scripts/hpc_cpu_matrix.py --total-cores 4 --ranks 1 2 4` 使用可選
OpenMP CLI，比較 `1×4`、`2×2`、`4×1`。每個配置首次執行另列，至少
三次正式重複，配置輪替，保留完整速度／壓力比較、geometry preflight、
rank RSS、phase 與啟動時間。實際核心綁定從每 rank 的 affinity 與 Linux
CPU topology 核對，SMT siblings 不算不同實體核心；亦要求原生組裝報告
顯示對應 team／batch。未量測模式不能因 build flags 而宣稱並行。

相關 Python 測試共 50 項通過，含錯誤核心預算、重疊綁定、缺少 worker
報告、錯誤 team 與 timeout 的拒絕。

`outputs/hpc02/hybrid/isolated-matrix/` 已完成 12 次執行；所有正確性
作業結束並確認資源空閒後才開始，期間未啟動其他編譯或模擬。每個配置
排除首次執行，以下為三次正式重複的中位數。所有 rank 的 phase、affinity
與原生 team 報告重新核對；三配置確實使用相同的 physical cores 0–3。

| ranks × threads | 端到端 (s) | 最大 rank 組裝 (s) | 最大 rank setup (s) | 最大 rank solve (s) | 最大 rank RSS (MiB) | 各 rank RSS 峰值總和 (MiB) |
|---|---:|---:|---:|---:|---:|---:|
| 1 × 4 | 14.451171 | 11.467306 | 2.286163 | 0.009694 | 253.582 | 253.582 |
| 2 × 2 | 20.822326 | 18.512145 | 1.745674 | 0.005742 | 193.074 | 378.617 |
| 4 × 1 | 34.692893 | 33.420429 | 0.757138 | 0.005140 | 143.824 | 556.418 |

這個案例的 `1×4` 與 `2×2` 相對純 MPI `4×1`，端到端分別快
2.400698 倍與 1.666139 倍。這是固定四核心的配置比較，不是跨節點 scaling
或所有案例的效能承諾。記憶體欄位是各程序各自峰值的總和，不是同時發生的
總峰值；也不能據此說較少 ranks 的每個程序都使用較少記憶體。

12 次原生收斂／質量 gate、geometry preflight 與配置間場比較均通過。
另外將全部輸出與 HPC-00C 封存的單 rank 參考比較，24 個速度／壓力比較
皆通過，最大相對 L2 為 `4.293576398446408e-15`。來源、輸入及 binary
在量測後重新核對並封存；`accepted-evidence.json` 綁定 summary、archive
及獨立參考比較。這份 accepted 僅指貼體固定四核心矩陣，FSI 正式量測
另外追蹤。

## 同 rank 單 thread 記憶體基準

最後稽核確認，固定總核心數比較同時改變 rank 數；例如 `1×4` 的最大
單 rank RSS 高於 `4×1`，即使各程序峰值總和較低，也不能直接宣稱所有
記憶體範圍符合改善門檻。為辨識 OpenMP 本身的增量，另以相同來源與輸入
重建無 OpenMP CLI，執行 1／2／4 ranks、各 rank 單 thread 的首次加三次
正式隔離基準，輸出至 `outputs/hpc02/hybrid/same-rank-serial-baseline/`。

比較方式預先固定：各混合配置對應相同 rank 數的單 thread 配置；同時報告
三次正式重複的最大「每 rank peak RSS 最大值」及最大「各 rank peak RSS
總和」之比，沿用 1.25 上限，保留總和不是同時總峰值的限制。固定四核心
的配置比較仍保留，不能用此基準抹除較少 ranks 的單一程序記憶體較高之事實。
數值仍須對單 rank 參考通過既定速度／壓力與原生守恆門檻。

同 rank 基準的 12 次執行已完成。其與固定四核心矩陣共 24 次執行的
原生 profile、配置、來源雜湊及 48 個封存單 rank 場比較已重新核對。
以下為三次正式重複取最大值後的比值：

| ranks | 混合配置 threads／rank | 最大單 rank RSS 比值 | 各 rank peak RSS 總和比值 |
|---:|---:|---:|---:|
| 1 | 4 | 1.008888 | 1.008888 |
| 2 | 2 | 1.005751 | 1.004995 |
| 4 | 1 | 1.001086 | 0.997138 |

六個比值都通過 1.25 預算。完整數值與原始 bytes 在
`same-rank-serial-baseline/memory-comparison.json`；來源、binary、建置命令、
輸入與量測證據已封存，`accepted-evidence.json` 記錄檔案雜湊。
核對程式為 `outputs/hpc02/hybrid/verify_same_rank_memory.py`。
結合 [FSI 效能及無 OpenMP 回歸](HPC_02_FSI_PERFORMANCE.md)，HPC-02 的
必要驗收已完成。這不表示後續分散式 FSI 或跨節點階段已完成。
