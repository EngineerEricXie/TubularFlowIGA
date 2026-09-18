# Latest goal objective / user override

請先閱讀 docs/plans/BIFURCATION_FSI_GPT_SOL_HANDOFF.md，再完整閱讀它連結的
PSC 時間效率改善計畫，並實際開始執行。从 P0 開始，依驗證門檻逐步完成；
保留現有未提交修改。所有計算與大型編譯透過 PSC Slurm 執行。
請實際進行修改與驗證，不要只重新提出計畫。

提交 Slurm 工作後，若沒有其他可執行的工作，使用**此對話的 sleep** 等待，
而不是背景監測腳本或結束本輪等待外部事件。Sleep 時長不固定為五分鐘：
根據已有耗時、剩餘 workload、排隊／執行狀態與 Slurm walltime 評估，
可以向上調整到二十分鐘、半小時、一小時或更久，以減少模型喚醒和 token 消耗。
每次 sleep 返回後，只檢查特定 Job ID 的必要 scheduler 狀態；
執行中不重複讀取相同求解日誌，不為「仍在執行」而做額外推理／回報。
只有作業完成、失敗、取消或超時後，讀取新 evidence 並自主接續下一張任務卡。
Sleep 被新使用者訊息中斷時，以最新指示為準。

完整目標與主計畫的驗證／交付要求不變。未通過 gate 的候選不得升級，
不可因等待或局部通過而將全目標標成完成。

此檔是使用者最新需求的持續交接記錄，並非宣稱已改寫產品內部 goal 主文字。
本會話只有 create_goal/get_goal/update_goal(status) 工具，沒有 objective 編輯介面；
不能用假完成／重建 goal 來繞過這項限制。
