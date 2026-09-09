# HPC-00C：單 GPU 與單程序基準

本報告接續 [CPU MPI 矩陣](HPC_00C_PROGRESS.md)，記錄單 GPU、浸入式與 FSI
路徑的首次／重複量測。完整清單狀態仍以
[WORKSTATION_HPC_TODO.md](../WORKSTATION_HPC_TODO.md) 為準。

## 工具與驗收範圍

新增 [hpc_serial_matrix.py](../../scripts/hpc_serial_matrix.py)，沿用原有
rank 包裝器、profile 彙整及場比較。每個案例執行一個首次程序及三個新的
重複程序；首次保留但不納入統計，不將其稱為已清除快取的測試。
每次程序 pin 在 logical CPU 0，OpenMP 與 BLAS threads 固定為 1。
不同案例依序執行；編譯、輸入準備與獨立驗證不與求解量測重疊。

GPU 模式先核對已接受 CPU 矩陣的輸入、二進位、量測工具及單 rank 首次
參考場雜湊，再執行 `device-info` 與實際求解。每個輸出都與該 CPU 參考
比較，使用原有 CPU/CUDA `1e-5` 相對 L2 與零參考 `1e-12` 絕對門檻。
速度與壓力分別比較，不調整 gauge；傳輸檢查完整節點 ID。

浸入式執行九組 centered-FD block 檢查，並沿用原生幾何、Newton 與守恆
門檻。FSI 沿用原生完整強耦合、合力／力矩投影、移動質量、壁面漏流與
離散連續性門檻。工具要求成功退出、完整單 rank profile 與有限診斷；
列印值有捨入，不取代原生程式中的完整精度判斷。

每次量測與中位數、min/max、母體標準差、各階段 exclusive 時間、host RSS
均保存。GPU 另要求唯一的 project-buffer peak 紀錄及結束時 live bytes 為 0。
失敗會保留已嘗試紀錄、停止矩陣，且不產生成功的重複統計。

## 建置、硬體與重現命令

執行環境延續 CPU 矩陣：WSL、Intel Core i9-14900KF，虛擬機可見 8 cores／
16 logical CPUs；GPU probe 確認 NVIDIA GeForce RTX 4080 SUPER，compute
capability 8.9、16 GiB device memory。CUDA 使用 `tubularflow-cuda` 環境與
SM89 目標。背景工作站及 Windows host 負載未受控，故不宣稱專用叢集效能。

```bash
make -C solvers/cpu immersed_aneurysm_jacobian_test compliant_channel_fsi_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
conda run -n tubularflow-cuda make cuda CUDA_ARCHS=89
python3 -m unittest discover -s scripts/tests -p 'test_hpc*.py' -v

conda run --no-capture-output -n tubularflow-cuda env \
  LD_LIBRARY_PATH=/home/tsungyeh/anaconda3/envs/tubularflow-cuda/targets/x86_64-linux/lib \
  python3 scripts/hpc_serial_matrix.py --case cuda-transport \
  --cpu-matrix outputs/hpc00/matrix/transport \
  --output-dir outputs/hpc00/serial-matrix/cuda-transport \
  --repetitions 3 --timeout 180 --cpu 0

conda run --no-capture-output -n tubularflow-cuda env \
  LD_LIBRARY_PATH=/home/tsungyeh/anaconda3/envs/tubularflow-cuda/targets/x86_64-linux/lib \
  python3 scripts/hpc_serial_matrix.py --case cuda-flow \
  --cpu-matrix outputs/hpc00/matrix/flow \
  --output-dir outputs/hpc00/serial-matrix/cuda-flow \
  --repetitions 3 --timeout 240 --cpu 0

python3 scripts/hpc_serial_matrix.py --case immersed \
  --case-dir outputs/hpc00/serial-matrix/immersed-inputs \
  --output-dir outputs/hpc00/serial-matrix/immersed \
  --repetitions 3 --timeout 300 --cpu 0

python3 scripts/hpc_serial_matrix.py --case fsi \
  --output-dir outputs/hpc00/serial-matrix/fsi \
  --repetitions 3 --timeout 1500 --cpu 0
```

重跑必須使用新的 output 目錄。浸入式輸入從先前已登記的 depth-2 副本
複製，來源 depth-3 fixture 保持不變。這些命令使用實際本機 GPU／MPI
socket 執行權限。此批重新編譯兩個 CPU fixture；CUDA 目標已是最新。
MPI 編譯器探測在 sandbox 內的 socket 訊息不是編譯器警告或數值失敗。

原始結果在 `outputs/hpc00/serial-matrix/`；各 `matrix.json` 保存當次
二進位／輸入／工具雜湊，`rank-0/run.json` 保存實際 affinity、命令及資源。
本批沒有改動數值核心、物理參數、檔案格式或原有驗收門檻。

## 單 GPU 已完成結果

兩種 GPU 案例共八次執行均接受。以下為子程序從啟動到結束的時間，包含
CUDA 初始化、求解與輸出；Python／taskset 外層啟動另有 launcher 計時。

| 案例 | 首次程序（s） | 三次重複中位數（s） | 重複範圍（s） | Host RSS 中位數（MiB） | Project device peak（bytes） |
|---|---:|---:|---|---:|---:|
| 流場 | 18.3552 | 17.0449 | 16.6405–17.0990 | 288.92 | 52,207,740 |
| 傳輸 | 0.8167 | 0.7178 | 0.7147–0.7645 | 282.57 | 24,744,428 |

GPU requested peak 在各次重複均相同，所有程序的 project-buffer live bytes
在結束時都是 0。這不包含 CUDA driver／library 的配置，不等於整張卡的
使用峰值；也不能把不同時刻的 host RSS 與 device peak 當成同步總記憶體。

| 案例 | 組裝中位數（s） | 預條件器建立（s） | 線性求解（s） | 輸出（s） |
|---|---:|---:|---:|---:|
| 流場 | 3.602011 | 0.001103 | 12.248704 | 0.003892 |
| 傳輸 | 0.045541 | 0.000048 | 0.359217 | 0.006414 |

表格使用 exclusive 時間，完整 phase 與每次樣本均在 summary 中。GPU
流場仍以線性求解為主要成本；原有 CPU 小案例則由組裝支配。兩後端採用
各自既有求解器，時間觀察不代表已完成預條件器擴展性工作。

八次 GPU 輸出的最大 CPU 相對 L2：流場速度 `4.80008e-8`、壓力
`2.74249e-8`、傳輸 `4.06867e-6`，均通過 `1e-5`。
矩陣結束後另對四份 GPU 流場執行 `iga_flow_validate`，檢查有限值及相對
質量不平衡 `<=1e-6`；最大值 `3.43214e-8`，四份全部通過。
命令、數值與日誌雜湊在 `independent-cuda-flow.json`。

CPU MPI 矩陣使用 core mask（可能包含兩個 SMT siblings），這批 serial/GPU
程序只 pin 一個 logical CPU；外層啟動路徑也不同。比較時應選擇一致的
子程序／application 計時範圍並註明配置差異，不直接拿不同 launcher
範圍的數字作速度比。GPU 傳輸 host RSS 高於 CPU 小案例，不能只報時間改善。

## 浸入式已完成結果

Depth-2 完整基準的首次與三次重複均通過。首次子程序約 `64.781 s`；
三次重複為 `65.0835、65.0870、65.5417 s`，中位數 `65.0870 s`，母體標準差
`0.2152 s`。Host RSS 中位數 `61,849,600 bytes`（約 `58.98 MiB`）。

Exclusive 組裝中位數 `62.949719 s`、幾何 `1.365466 s`、預條件器建立
`0.179712 s`、線性求解 `0.001980 s`。這是含九組 FD block 檢查的完整
回歸時間，不是單次求解時間；四次都保留原生完整 gate。

最大列印 FD defect 是 `1.41286e-11`，最終 Newton residual
`2.50713e-17`，open-port normalized balance `2.08177e-15`，wall leakage
normalized `1.573e-5`。各次完整數值診斷在 `immersed/summary.json`。
此單程序結果不構成浸入式 MPI 或 OpenMP 支援證據。

## FSI 已完成結果

完整首次與三次重複均正常退出並通過原生驗收。首次子程序 `634.168858 s`；
重複樣本為 `632.715026、629.447414、630.371758 s`，中位數 `630.371758 s`，
母體標準差 `1.375282 s`。Host RSS 中位數 `75,399,168 bytes`（約 `71.91 MiB`）。

Exclusive 組裝中位數 `620.625380 s`、幾何 `7.910823 s`、預條件器建立
`0.684716 s`、線性求解 `0.007759 s`、耦合自身工作 `0.721206 s`。
耦合 inclusive 時間包含其內部組裝等階段，不可再次與這些子階段相加。
這些數字涵蓋完整基準的組裝呼叫，不只最後一個 Newton 步。

四次皆以四個強耦合迭代完成；最終列印位移殘差 RMS `2.71104e-8 m`，
門檻 `4.80669e-8 m`。中心位移 `3.80669e-5 m`；moving mass 與 wall leakage
皆為 `0.0155613`，discrete continuity `1.82225e-14`。原生完整精度的
合力／力矩各分量 `1e-11` 限制及其餘檢查均執行通過。這仍是單分區 FSI。

## 測試與完成範圍

30 項 Python 測試通過。新增測試拒絕缺少／重複／非有限 FD、缺失的 FSI
收斂歷史、錯誤 CUDA allocation scope、未釋放的 project buffer、變更過的
CPU 參考場／輸入，以及失敗／缺少／重複的量測樣本；首次不納入重複統計。

另將獨立測試副本的 `surface.vtp` 截斷，執行相同 immersed 矩陣工具。
原生 fixture 失敗後 controller 退出 1，只保留一次失敗的首次程序，未執行
重複樣本，`repeated_statistics` 為 null；紀錄位於 `bad-immersed/`。
十二份成功 profile 與三組重複統計均另行重算核對；GPU 流場的四次獨立
質量檢查也全部通過。

FSI 結束後，十六份 serial／GPU profile、四組統計、實際 CPU affinity 及
輸入／二進位／量測工具雜湊全部重新核對；另核對先前 CPU／serial 證據的
1,061 個檔案雜湊及四份獨立 GPU 質量驗收。最終彙整在
`outputs/hpc00/serial-matrix/accepted-evidence-complete.json`，含本批 186 個
產物雜湊。較早的 partial 報告與失敗樣本原樣保留。

HPC-00C 的選定單機模式矩陣已完成：CPU MPI 32 次，加上 serial／GPU
16 次。尚未實作的 OpenMP／混合及分散式 immersed／FSI 模式為 N/A。
量測到的記憶體增加與負收益仍保留，不宣稱所有配置均達成效能改善目標。
HPC-00D 的獨立物種收支、HPC-02／03／07 的平行實作及 HPC-09 的跨節點
驗收仍依各自任務追蹤。
