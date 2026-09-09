# HPC-00D：數值參考與效能門檻

狀態：選定四個基準的驗收方式、數值參考與優化門檻已建立並驗證，
HPC-00D 完成。此結論不代表後續 OpenMP／MPI／分散式 FSI 已完成。
執行日期：2026-09-08 UTC；工作區包含未提交修改，參考依照記錄的
source／binary／input 雜湊識別，不能僅以 Git HEAD 識別。

## 固定數值門檻與權威來源

通用比較政策位於 [hpc_baselines.json](../../benchmarks/hpc_baselines.json)，
物種收支與完整 fixture 狀態分別補充於
[hpc_transport_budget.json](../../benchmarks/hpc_transport_budget.json) 及
[hpc_reference_states.json](../../benchmarks/hpc_reference_states.json)。
三份政策都在各自驗收與後續效能優化之前登記；沒有放寬失敗的門檻。

| 案例 | 數值參考與全場比較 | 獨立或原生必要 gate |
|---|---|---|
| 貼體流場 | 相同生成輸入的 CPU 單 rank MUMPS 參考；速度與壓力分別比較；CPU `1e-6`、CUDA `1e-5` 相對 L2 | 正 Jacobian、Newton 殘差、有限值、相對質量不平衡 `<=1e-6` |
| 貼體物種傳輸 | 相同生成輸入的 CPU 單 rank 原生 GMRES 參考；初始及兩個時間步；CPU `1e-6`、CUDA `1e-5` 相對 L2 | 正 Jacobian、原生 KSP 成功、各物種體積／表面收支、自由列殘差及指定濃度 `1e-6`；淨線性轉換零參考 `1e-12` |
| 浸入式 depth 2 | 完整原生 fixture 的 216 節點速度／壓力與一個 controller；逐物理量 `1e-6` 相對 L2 | 九組 centered FD `<=1e-8`，零 action 絕對門檻 `1e-11`；Newton `<=1e-11`、首步真實線性殘差 `<=1e-10`、open／wall 守恆 `<=1e-3`；原生 cap 面積／法向等幾何 gate |
| compliant-channel FSI | 完整原生 accepted macro-step，九個物理量各自比較 `1e-6` 相對 L2 | 原生位移 RMS 強耦合門檻、移動幾何與提交一致性、moving mass／wall leakage `<0.03`、continuity `<1e-8`；合力差 `<=1e-11 N`、力矩差 `<=1e-11 N m` 各分量 |

所有全場比較都使用明確節點 ID 與欄位定義；零參考改用絕對 L2 `1e-12`。
不移除 pressure offset、不調整物理單位。這些是係數範數；物種收支與
FSI 原生面積加權 RMS 是另外的物理／離散方程 gate，不能互相替代。
幾何驗收限於所選 fixture 與實際積分／取樣範圍，不宣稱所有可能幾何都有效。

FSI 力／力矩完整精度的權威判斷來自原生測試中的 `Require`；輸出摘要
存在捨入，不能據摘要重建完整精度 gate。完整狀態輸出在所有原生 gate
通過後才發布。參考是數值回歸基線，不是實驗資料或臨床有效性證明。

## 已核對證據

- [CPU MPI 矩陣](HPC_00C_PROGRESS.md)：1／2／4／8 ranks，流場及傳輸共
  32 次；8 份 packed-geometry 檢查及 16 份獨立流場質量檢查。
- [單 GPU 與 serial 矩陣](HPC_00C_SERIAL_PROGRESS.md)：16 次選定模式，
  四份獨立 GPU 流場質量檢查；FSI／浸入式原生 gate 全部執行。
- [物種收支](HPC_00D_TRANSPORT_BUDGET.md)：69 項解析／錯誤檢查、
  CPU 1／2／4／8 rank 與 CUDA 共 10 個時間步，以及五種 CLI 錯誤驗收。
- [完整浸入式／FSI 狀態](HPC_00D_REFERENCE_STATES.md)：四個獨立原生
  collection、12 個物理量比較全部相對 L2 為 0；22 項 writer 檢查、
  41 項 Python 測試、七種原生拒絕及實際 CLI 數值擾動拒絕。

每份報告含重現命令、來源、硬體／軟體、執行方式、原始檔案位置與限制。
CPU 參考由相同輸入重新產生，原生 fixture 由其完整定義及原生 gate
確定。最新完整 fixture 來源與執行檔另封存供追溯；後續改版比較使用
已接受的舊參考，不要求舊參考的 binary 路徑仍指向當時執行檔。
新的執行仍必須記錄自己的 source／binary／inputs 並核對案例身分。

機器可讀證據在 `outputs/hpc00/matrix/`、`serial-matrix/`、
`transport-budget/` 與 `reference-states/`。各目錄的 accepted evidence／
inventory／matrix／collection 記錄組成參考來源鏈；生成檔案不提交。
較早的 partial 報告與失敗觀察保留為歷史，不能將其狀態誤讀為目前結論。

## 效能目標與未完成範圍

後續優化的既定改善宣稱要求：至少三個獨立新程序重複，以端到端中位數
改善至少 10%，同時數值 gate 不退步。RSS 比例門檻是 `1.25`，必須同時
說明最大 rank 與總和範圍，以及不同時間峰值相加的限制。
組裝、求解、預條件器建立、通訊與輸出分別記錄；GPU peak 另有範圍定義。

目前小型 MPI 案例已量測到記憶體增長超過門檻，報告保留負收益，沒有
將基準建立解讀成所有效能條件通過。預處理成本、OpenMP 多核心組裝、
浸入式／FSI 分散計算、耦合續跑、大型 I/O 與跨節點 scaling 繼續由
HPC-01 至 HPC-09 各項驗收；本報告不縮減那些工作的完成要求。
