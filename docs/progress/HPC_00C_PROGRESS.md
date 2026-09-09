# HPC-00C：單機重複執行矩陣進度

本批補上貼體 CPU 流場／傳輸的純 MPI 矩陣工具與實測。後續單 GPU 與
浸入式與 FSI 量測見 [單程序進度](HPC_00C_SERIAL_PROGRESS.md)。選定模式的
完整首次／三次重複矩陣已完成，HPC-00C 已勾選；不代表後續平行化已完成。
OpenMP／混合組裝目前未實作，應標為 N/A，不以增加環境變數宣稱支援。

## 交付與執行方式

新增 [hpc_cpu_matrix.py](../../scripts/hpc_cpu_matrix.py)，沿用每 rank 包裝器、
profile 彙整及場比較工具。每個配置先執行獨立幾何檢查，再執行一次首次
求解及三次 fresh-process 重複。重複輪次輪換配置順序，各求解依序執行。
首次求解不納入時間統計，也不代表檔案快取已清除。

工具保留完整 rank 記錄、數值比較、時間樣本、中位數、範圍、母體標準差、
階段時間與 RSS。輸入、二進位及量測腳本在執行前後核對雜湊。失敗會保留
所有已嘗試紀錄並停止；不刪除失敗樣本後繼續宣稱矩陣通過。

本批沒有修改 C++ 求解器、物理參數、檔案格式或數值容許值。
基準為現有工作目錄版本，不能只用 HEAD 表示已有大量未提交修改的程式。
生成資料及完整結果保留於 `outputs/hpc00/matrix/`；不提交這些產物。

```bash
make -C solvers/cpu iga_solve iga_navier_stokes iga_mesh_check iga_pack \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
python3 -m unittest discover -s scripts/tests -p 'test_hpc*.py' -v

python3 scripts/hpc_cpu_matrix.py \
  --case-dir outputs/hpc00/matrix/transport-inputs \
  --database-stem straight_neurite --system neuron_transport \
  --output-dir outputs/hpc00/matrix/transport \
  --ranks 1 2 4 8 --repetitions 3 --timeout 180

python3 scripts/hpc_cpu_matrix.py \
  --case-dir outputs/hpc00/matrix/flow-inputs \
  --database-stem straight_tube --system flow \
  --output-dir outputs/hpc00/matrix/flow \
  --ranks 1 2 4 8 --repetitions 3 --timeout 240
```

重跑須換新的 output 目錄。案例沿用 HPC-00A/B 的 1005-node、720-element
生成資料，複製至獨立目錄後，以 `mpmetis` 產生 4／8 分區，再用 `iga_pack`
打包。既有 1／2-rank 資料庫保留。準備命令及退出碼見 `preparation.json`，
每個配置的幾何檢查與求解分開計時。

本機是 WSL 上的 Intel Core i9-14900KF，虛擬機 topology 顯示 8 cores、16 logical CPUs；Open MPI 使用
`--bind-to core --map-by core --nooversubscribe`。每 rank 的實際 affinity
另存在 `rank-N/run.json`，每個 core 的兩個 logical CPU 均可能出現在 mask 中。
OpenMP 與 BLAS threads 都固定為 1。流場使用已登記的 MUMPS 直接參考；
傳輸保留既有預設求解器。MPI 實測使用允許本機 socket 的執行權限。
本批沒有跨節點排程作業，因此不構成 HPC-09 的驗收。
量測期間沒有同時啟動其他模擬或編譯；工作站仍有背景程序，Windows host
負載與檔案快取未受控。以下是此環境的觀察值，不是專用節點的效能承諾。

## 傳輸結果

四個配置的幾何檢查及全部 16 次求解／場比較通過。以下中位數只使用三次
重複執行；端到端時間包含 MPI／Python 啟動、求解、輸出與結束。

| Ranks | 端到端中位數（s） | 相對單 rank 速度比 | 最大 rank RSS 中位數（MiB） | RSS 峰值總和／單 rank |
|---:|---:|---:|---:|---:|
| 1 | 1.8663 | 1.000 | 89.55 | 1.000 |
| 2 | 1.7159 | 1.088 | 75.01 | 1.578 |
| 4 | 1.6164 | 1.155 | 66.12 | 2.728 |
| 8 | 1.6166 | 1.154 | 60.23 | 4.919 |

每個 rank 峰值的總和不是同一時刻的總 RSS。多 rank 的此項比例超過
預先登記的 1.25 記憶體門檻，不宣稱所有效能條件通過。4／8 rank 在這個
小案例的時間相近；單次啟動成本及量測包裝器也占端到端時間，不能外推大型
案例或跨節點效率。完整變動範圍與各次樣本在 `transport/summary.json`。

相對新產生的單 rank 首次參考場，2／4／8 rank 的最大相對 L2 分別為
`1.81155e-7`、`1.69284e-7`、`3.69511e-7`，均低於不變的 `1e-6` 門檻。
比較檢查完整場形狀、節點 ID 與有限值。原生求解成功與場比較不取代
獨立物種守恆驗收；後續結果見 [物種收支報告](HPC_00D_TRANSPORT_BUDGET.md)。

| Ranks | 組裝（s） | 預條件器建立（s） | 線性求解（s） | 顯式通訊區段（s） |
|---:|---:|---:|---:|---:|
| 1 | 1.301981 | 0.023992 | 0.014656 | 0.000702 |
| 2 | 1.169276 | 0.004634 | 0.021166 | 0.029001 |
| 4 | 1.132612 | 0.000752 | 0.014483 | 0.346665 |
| 8 | 1.058148 | 0.000322 | 0.008572 | 0.563236 |

表格是每次最大 rank exclusive 時間的中位數；各區段的最大值可能來自
不同 rank，不能相加成 critical-path 時間。PETSc 內部通訊仍屬呼叫所在區段。

## 流場結果

四個配置的幾何檢查及全部 16 次求解／速度／壓力比較通過。每次首次與
重複執行均保留在 `flow/summary.json`；統計仍只取三次重複。

| Ranks | 端到端中位數（s） | 重複範圍（s） | 相對單 rank 速度比 | 最大 rank RSS（MiB） | RSS 峰值總和／單 rank |
|---:|---:|---|---:|---:|---:|
| 1 | 35.2654 | 35.2644–35.3621 | 1.000 | 251.56 | 1.000 |
| 2 | 32.0563 | 31.7522–32.3589 | 1.100 | 192.36 | 1.497 |
| 4 | 31.0558 | 31.0065–31.3607 | 1.136 | 143.91 | 2.203 |
| 8 | 30.0182 | 29.4635–30.3142 | 1.175 | 118.58 | 3.384 |

速度場最大相對 L2 是 `4.29358e-15`，壓力是 `1.34037e-15`；均低於既定
`1e-6`，沒有調整壓力 gauge。所有求解以兩次 Newton 更新收斂；單 rank
初始殘差 `0.0735831`、最終殘差約 `8.4761e-15`，收斂容許值
`7.35831e-10`。精確逐次診斷保留於 rank 日誌。

| Ranks | 組裝（s） | 預條件器建立（s） | 線性求解（s） | 顯式通訊區段（s） |
|---:|---:|---:|---:|---:|
| 1 | 32.559928 | 2.043140 | 0.009011 | 0.001381 |
| 2 | 29.999724 | 1.536871 | 0.005767 | 0.387923 |
| 4 | 29.896688 | 0.684099 | 0.004386 | 8.834329 |
| 8 | 28.977293 | 0.601272 | 0.004784 | 15.386060 |

通訊區段包含等待；例如較早完成組裝的 rank 可以在其他 rank 仍組裝時等待。
因此 8-rank 的組裝最大值與通訊最大值不能相加。整體仍受組裝支配。
8-rank 端到端時間減少約 14.9%，但 RSS 峰值總和增加至 3.38 倍，超出
1.25 門檻；兩個案例都不能宣稱時間與記憶體條件全面通過。

矩陣結束後，另對全部 16 份流場輸出執行
`iga_flow_validate DATABASE FIELD`，檢查回傳碼、所有輸出量的有限性與
不變的 `1e-6` 相對質量不平衡門檻。最大質量不平衡為 `3.43214e-8`，全部
通過；最大相對散度定理誤差 `1.67600e-8` 另作診斷，未冒充新增的驗收門檻。
命令、數值與日誌雜湊保留在 `independent-flow-summary.json`。八次幾何
preflight 均為 720 元素、最小 detJ `9.01032e-5`、零 bad elements/samples。

## 分區工作量的證據與後續方向

`partition-load.json` 讀取本次 v5 資料庫的 element owner 與 required-element
index；兩個案例使用相同幾何，得到相同分布：

| Ranks | 每 rank required 元素數 | 全部 required 次數／720 |
|---:|---|---:|
| 1 | 720 | 1.000 |
| 2 | 658、673 | 1.849 |
| 4 | 469、658、673、458 | 3.136 |
| 8 | 308、451、534、580、631、513、452、296 | 5.229 |

8-rank METIS owner 元素數是 `87、92、92、92、92、91、87、87`。
owner 分配接近平均，但它不等於矩陣組裝工作量。
[Database::NodeRange](../../solvers/cpu/include/IgaDatabase.hpp) 依節點編號
連續分配矩陣列；[OwnedRowAssembler](../../solvers/cpu/include/OwnedRowAssembler.hpp)
載入 touching/required 元素。各 rank 計算局部元素矩陣後，只插入自己擁有
的列。因此同一元素可在多個 rank 積分，而全域矩陣列不會因此重複累加。

這些資料支持優先量測與改善 required 元素重複積分及資料配置的方向。
後續 HPC-06A／HPC-02 應比較內部節點重排、計算元素所屬 rank 後交換貢獻、
或保留 owned-row 組裝並平行化局部積分的成本。任何重排須保留外部節點 ID、
速度／壓力輸出、checkpoint 與 `.ntiga` 相容性；不能僅更改 owner 欄位。
本批只有量測及分析，尚未實作這些優化。

## 測試與剩餘工作

24 項 Python 測試通過，涵蓋首次樣本排除、失敗／缺少／重複樣本拒絕、
不一致 phase schema、非有限量測、timeout、非零退出、輸出目錄不可覆寫，
以及既有來源追蹤、每 rank 包裝器與場比較測試。受影響的 CPU 目標由 Make
確認已是最新；實際 MPI 數值結果另以上述執行紀錄為準。

另在獨立副本中故意以一分區資料庫取代兩分區檔案，執行相同矩陣工具。
一 rank 幾何檢查通過、兩 rank preflight 有限時間失敗；controller 退出 1，
`wrong-partition/summary.json` 保留兩次 preflight、零次求解、失敗原因，
`repeated_statistics` 為 null。全部 32 份成功 profile 與重複統計另重新
計算核對通過，沒有用缺少或失敗樣本作出成功統計。

後續已補上浸入式與單 GPU 的完整重複量測，FSI 尚待完成；見
[單 GPU 與單程序基準進度](HPC_00C_SERIAL_PROGRESS.md)。
HPC-00D 的獨立物理門檻與 HPC-09 的真實跨節點證據也保持各自待辦。
