# gpt-sol 執行交接入口

本輪只完成計畫，尚未執行新一輪優化。收到使用者的執行指示後，從 P0 開始；不要直接再跑一次 20 步心跳。

完整計畫：[PSC 時間效率改善計畫](BIFURCATION_FSI_PSC_TIME_OPTIMIZATION_PLAN.md)。

## 可直接貼給下一個 session 的指令

> 請執行 `docs/plans/BIFURCATION_FSI_PSC_TIME_OPTIMIZATION_PLAN.md`。目標是在相同問題與驗證要求下，縮短 PSC 上取得完整 FSI 結果的時間，並報告 CPU-hours。先完整讀計畫及最新專案指示，核對現有 jobs 和未提交修改。第一輪完成 P0 的來源凍結、可重跑測試入口及細項計時，再處理 P1-A 的 FSI 前置重複組裝；每張卡完成後依 gate 自主接續。保留現有工作，不重做已通過的平行組裝、VTK 修復或舊完整模擬。計算與大型編譯透過 Slurm，在 node-local scratch 執行，保存 immutable evidence。每張卡更新 STATE.md/state.json；未通過驗證的候選不能用在正式長跑。遇到無收益分支依計畫跳過並記錄理由，不無限測試。提交或發布依當時授權。

## 首輪具體工作

1. CWD：`/ocean/projects/mch260002p/thsieh1/TubularFlowIGA`。
2. 讀完整計畫第 1、2、3、6、7 節；確認當前 AGENTS.md／使用者指示。
3. 查 `git status --short`、`git diff --stat`、`squeue -u "$USER"`。不自動取消其他作業。
4. 核對 `benchmarks/heartbeat-optimized-45982619/` 與 `benchmarks/heartbeat_assembly_comparison.json`。
5. 建立 `benchmarks/fsi-time-opt/<run-id>/STATE.md`、`state.json`；先標 P0-A running，其餘 pending。
6. 凍結目前來源：包含既有 dirty/untracked source，排除大型結果／環境／憑證。不只封存 HEAD。
7. 確認 archived tested source 和目前 source 的數值等價性；目前工作樹不能不加檢查就當成可比 B0。
8. 在 Slurm 內整理 build/test runner，W0/W1 小測可重現之後，再加呼叫原因／細項計時。
9. P0 gate 通過後實作 P1-A；先證明 `BeginTrial→Assemble→SolveTrial` 中的第一個 Assemble 沒有必要的副作用，再移除，跑目標 regression。
10. 每個階段都保留 source/binary/input hash、Slurm JobID、原始 logs、正確性和時間數據。下一張卡看主計畫，不憑這份摘要跳過測試。

## 最容易踩錯的地方

- `45982619`：20 步求解成功，耗時 13603 秒；job FAILED 是最後 VTK metadata 驗證失敗。
- `46006012`：已修好 metadata、native VTK 驗證和打包成功；不用重新修一次。
- `heartbeat-paraview.tar.gz` 已存在；基準仍是 `visualization_only`，17 步未過原守恆門檻。
- 8 threads 固定幾何 assembly 比新版 1 thread 快 3.64 倍；不是全程 speedup，也不是 CPU 利用率。
- 現有 volume/trace/wall worker 平行化已通過，不要再當新方案實作。
- `IGA_ASSEMBLY_BATCH_SIZE` **已存在**，預設 `min(threads,8)`；大於 8 threads 時尤其要測 batch。
- `IGA_PROFILE_DETAIL`、`AssemblyRequest`、`run_phase.sh` 是**計畫提出、尚未實作**的名稱。
- 目前 heartbeat runtime 使用 SeqAIJ/COMM_SELF；增加 `mpiexec -np` 不是可用的加速方式。
- OpenMP 測試要 MPI_THREAD_FUNNELED，不能放寬 thread-safety gate；檢查既有 wrapper。
- VTU 不是完整 restart，不能直接從展示檔續算峰值。
- PSC RM-shared 不使用 `-C PERF`；用內部計時先定位，必要硬體 profiling 才另安排短 RM 作業。
- 沒有證據前不改 dt、積分深度、網格、材料、coupling／守恆門檻來取得加速。

首輪完成時回報：P0 evidence 路徑、call-count hot paths、P1-A 改动與測試、候選是否已可升級、下一個具體任務。不得只回覆「計畫可行」。
