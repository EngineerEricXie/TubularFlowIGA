# HPC-00D：獨立物種收支驗收

狀態：本報告的獨立物種收支驗收已通過。後續浸入式／FSI 完整狀態
參考亦已通過，整項完成範圍見 [HPC-00D 報告](HPC_00D_REPORT.md)。
執行日期：2026-09-08 UTC。HEAD 為 `ee8da2ab528849814b6d12c8186476b5f7f59124`，
另有既有及本批未提交修改；完整 source／binary 雜湊保存在下述產物。
此報告接續 [原有數值基準](HPC_00BD_PROGRESS.md)。

## 驗收範圍與方程

新增 `iga_transport_validate`，直接讀取 `.ntiga`、案例設定、指定速度場及
相鄰時間步的節點係數。每個資料庫元素恰好積分一次，不讀取求解器的
PETSc 矩陣、元素矩陣或殘差；MPI 各分區輸出的全域場使用相同程序檢查。
此離線工具目前將完整場與自由列殘差存於單程序記憶體，適用於選定的小型
正確性基準；大型分散式診斷的記憶體需求仍由 HPC-06 追蹤。

對每個物種，體積收支殘差為：

```text
R_volume = capacity_rate + advection_volume + reaction - source
           + natural_outward - essential_inward
```

`capacity_rate` 由 backward Euler 的前後濃度差與容量係數積分。
擴散使用弱式梯度項；指定通量與 Robin 條件計入 `natural_outward`。
原生 Flux 值在方程右側為正，因此工具將其取負作為向外通量。
Robin 向外交換量為 `coefficient * (current - exterior)`。
Dirichlet 的 `essential_inward` 則是施加約束前的弱式列反作用量，
不將它等同於直接由邊界濃度梯度計算的通量。

原生平流項使用 `v · grad(c)`，故以表面通量表達時必須保留速度散度修正：

```text
R_surface = capacity_rate + advection_outward - integral(c * div(v))
            + reaction - source + natural_outward - essential_inward
```

只有總量閉合仍可能漏掉局部相抵的錯誤，因此另檢查未受 Dirichlet 約束列的
殘差 L2，以及前後兩步的指定濃度。兩物種基準另要求淨線性轉換為零。
濃度積分的單位是濃度乘物理體積；沒有單位轉換依據時，不稱為公斤。

數值門檻於執行前登記於
[hpc_transport_budget.json](../../benchmarks/hpc_transport_budget.json)：
相對門檻 `1e-6`，零參考絕對門檻 `1e-12`。容量正規化保留前後兩個
inventory/dt，邊界使用未相抵的通量大小；各列殘差以各組成項 L2 之和
正規化。完整定義保存在政策檔，不以本次求解結果調整門檻。

支援固定貼體幾何、指定節點速度、常數係數、backward Euler 與純量場。
非擴散方程可使用既有 SUPG 測試函數。原生 SUPG diffusion 的離散梯度項
含額外測試函數乘因子，此工具明確拒絕該模式；所選 N0／Nplus 基準中，
N0 有擴散但沒有 SUPG，Nplus 有 SUPG 但沒有擴散，因此在支援範圍內。
動態邊界、移動幾何與耦合速度並非本工具目前的驗收範圍。

工具共用幾何／基底與設定解析元件，因此其獨立性限於場插值及方程收支
的重新積分；不能取代解析幾何、網格收斂或所有物理模型驗證。

## 已執行命令與證據

```bash
make -C solvers/cpu transport-budget-test

conda run --no-capture-output -n tubularflow-cuda env \
  LD_LIBRARY_PATH=/home/tsungyeh/anaconda3/envs/tubularflow-cuda/targets/x86_64-linux/lib \
  python3 scripts/hpc_transport_budget_regression.py \
  --cpu-matrix outputs/hpc00/matrix/transport \
  --output-dir outputs/hpc00/transport-budget/validation1 \
  --ranks 1 2 4 8 --include-cuda
```

輸出目錄必須尚未存在；CPU 參考矩陣的案例、二進位與量測來源雜湊必須
仍匹配。這批命令增加逐步輸出，不能與先前效能矩陣直接比較時間。

解析測試涵蓋體積尺度、非零速度散度、擴散通量方向、Robin 外部濃度、
跨物種轉換，以及全域積分相抵但自由列殘差不為零的錯誤。
原生回歸產生 CPU 1／2／4／8 rank 及 CUDA 的初始場與兩個時間步；
每一步均通過既有場比較門檻及逐物種收支與殘差檢查。
另外以自由係數擾動、非整數節點 ID、截斷檔案、非有限值及額外欄位
確認實際 CLI 以非零退出碼拒絕錯誤。

CLI 退出碼：`0` 通過所有通用數值門檻，`2` 數值門檻失敗，`1` 輸入或
執行失敗。淨物種轉換為零是此雙物種案例的附加驗收，由回歸驅動程式
強制檢查，不將它套用到所有可能含生成或衰減的泛用系統。

建置使用 GCC 11.4、C++17、`-O3 -Wall -Wextra -Wpedantic`，沒有新增警告；
解析／拒絕測試共 69 項通過，日誌為 `outputs/hpc00/transport-budget/build-unit.log`。
CPU 執行使用既有 PETSc 3.15.5／OpenMPI，單 thread、core binding，保持
transport 原生 GMRES 預設；CUDA 使用 RTX 4080 SUPER、SM89、Conda CUDA 12.6。
FSI 重複量測結束後才啟動此批建置及驗收。

| 模式 | 最大體積收支相對缺陷 | 最大表面收支相對缺陷 | 最大自由列相對殘差 | 最大場相對 L2 |
|---|---:|---:|---:|---:|
| CPU 1 rank | 8.96456e-8 | 8.96456e-8 | 2.90806e-8 | 0（參考場） |
| CPU 2 ranks | 9.37344e-8 | 9.37344e-8 | 7.19889e-8 | 1.81155e-7 |
| CPU 4 ranks | 9.13514e-8 | 9.13514e-8 | 4.74988e-8 | 1.69284e-7 |
| CPU 8 ranks | 9.81195e-8 | 9.81195e-8 | 6.51320e-8 | 3.69511e-7 |
| CUDA | 7.68025e-8 | 7.68025e-8 | 2.00596e-7 | 4.06867e-6 |

各列取兩個時間步、兩個物種的最大絕對正規化缺陷。全部十個時間步的
淨線性物種轉換為零，前後指定濃度檢查亦通過。每次獨立工具恰好積分
720 個元素；每份場包含 1,005 個節點。CPU 場門檻 `1e-6`、CPU/CUDA
場門檻 `1e-5` 及原登記收支門檻均保持不變。

自由濃度係數擾動退出 `2` 且自由列殘差 gate 失敗；非整數 ID、截斷、
非有限值及額外欄位均退出 `1`。回歸驅動程式整體退出 `0`，
`validation1/summary.json` 狀態為 `requested_checks_passed`。

執行後另核對輸入／binary／source 快照、五份 rank profile、所有逐步場比較
及十份原始收支 JSON，結果與摘要一致。`outputs/hpc00/transport-budget/accepted-evidence.json`
記錄最大缺陷及 188 個產物雜湊。原始日誌及機器可讀資料在 `validation1/`，
本批不修改生產數值核心或檔案格式。

此小案例結果不宣稱大型驗收工具的 MPI 擴展性，也不替代 FSI 分散式
實作或跨節點驗收。

另已核對兩個 CPU 矩陣的 32 份 profile／場比較、8 份正 Jacobian
檢查、16 份獨立流場質量檢查，以及 16 份 serial／GPU 的來源雜湊與原生
診斷，彙整於 `outputs/hpc00/baseline-acceptance.json`。該檔的
`selected_baseline_acceptance_verified` 僅表示已列明的基準檢查已核對，
不是整份 HPC 清單或整項 HPC-00D 完成標記。

本批之後已為浸入式與 FSI 補上完整且具有穩定 ID 的狀態參考，
登記相對 L2 比較方式並驗證拒絕錯誤的能力，見
[完整狀態參考](HPC_00D_REFERENCE_STATES.md)。現有原生殘差、守恆及
合力／力矩門檻繼續保留；成功退出與相同摘要值，不能代替完整場在後續
組裝平行化前後的比較。
