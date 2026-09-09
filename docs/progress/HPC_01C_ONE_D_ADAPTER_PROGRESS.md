# HPC-01C：1D adapter 的本地操作與求解準備

日期：2026-09-08。狀態：**本報告範圍通過；HPC-01C 整體仍部分完成。**

## 實作與契約

`OneDFlowDomainAdapter` 與 `OneDFlowTransportDomainAdapter` 現在透過 runtime
新增的 `RunLocalAdapterStage` 使用既有 `FailureAgreement`。BeginStep、port
輸入、PrepareCommitStep、Abort 與本地 rollback 都先捕捉例外，再交給程序群
協調。staged hydraulic／transport 準備及 rollback 前置驗證也納入協調，
避免某 rank 因缺少輸入離開，而其他 rank 已進入原生 MPI 求解或 rollback。

會使用 collective 的原生 Solve／staged rollback 保留在 local work callback
之外。PrepareCommitStep 全群組成功後才可呼叫 noexcept FinalizeCommitStep。
getter 與 constructor 保留本地語義；沒有為單一 rank 的 port 查詢加入 MPI。
建構與 executor 中 getter 後的本地處理，仍由 embedding 層負責錯誤協調。

staged `SetPortInput` 與 `SetTransportConcentration` 先建立候選 map，整個程序群
成功後才 swap。輸入被拒絕時不留下其他 rank 已保存的重複資料，修正後可重送。
transport 資料準備或 replay 失敗，依原有契約清除該 scalar attempt 的輸入，
重試時必須重新提供完整濃度；已完成的 hydraulic frames 可繼續使用。

transport 的 phase 驗證在 scalar attempt 外協調，保留原有語義：拒絕非法
重複求解不撤銷先前成功 trial，也不清除尚可使用的輸入。hydraulic 成功旗標
同樣只在準備成功、開始新求解時清除。測試驗證拒絕非法重複求解後，port 查詢
與提交仍可成功。

啟用 agreement 後，同一 group 必須以相同順序呼叫 adapter 的變更操作；
不可將這些已協調的方法再包進另一個 group-local work callback。未注入
agreement 的使用者保留純 C++ 行為。generic／bifurcation graph 的既有 runtime
已注入 communicator callback，因此這些 adapter 操作直接接上該協議。
其他自行建構 runtime 的 embedding 必須明確注入；本次不宣告所有舊 CLI 都已接入。

沒有更改物理公式、port 單位／方向、logical-to-native species binding、
coupling scheme、staged vasodilation 限制或檔案格式。候選 map 增加局部複製
成本，仍需在 HPC-00C／06A 的大型案例量測。

## 驗收

環境沿用 WSL、Intel i9-14900KF、GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5
real/double、32-bit PetscInt、MUMPS。HEAD
`ee8da2ab528849814b6d12c8186476b5f7f59124` 加未提交工作樹，精確來源見 inventory。
原生六物種 fixture 只在測試記憶體配置關閉 vasodilation，以符合既有 staged
能力限制；檔案不變。每次 trial 為三個 configured substeps。

```bash
make -C solvers/one_d core-test one_d_adapter_failure_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
make -C solvers/coupling petsc multidomain_subcommunicator_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
OMP_NUM_THREADS=4 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 \
  PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps' \
  timeout 180 mpiexec --oversubscribe -np 3 \
  solvers/one_d/one_d_adapter_failure_test examples/one_d/multispecies_physiology
```

同一 adapter test 另以 1 thread 執行。測試從實際 OpenMP region 計算 thread
數，三個 rank 的成功紀錄均符合要求。各 thread 配置都執行 world 及獨立
1+2 groups，分別通過 20、20、40 個故障案例；兩 rank group 的重複次數較多，
用來暴露意外的 world collective。數量以 group 計算，不將各 rank 重複計數。

20 種模式含 16 個本地階段的單 rank 拒絕，以及四個真實輸入錯誤：錯誤時間、
不完整 logical concentration、缺 hydraulic input、缺 inward concentration。

| 驗收 | 結果 |
|---|---|
| 本地故障協調 | 每個成員收到包含指定 rank 與階段的錯誤；進入求解前的故障皆沒有執行 implicit advance |
| 求解邊界 | 後段故障只完成原定三次 implicit advance，沒有進入額外求解 |
| staged 輸入重送 | mixed hydraulic／scalar map 被共同拒絕後可重送，不產生 duplicate-input 錯誤 |
| scalar 重試 | 修正或重送濃度可沿用已完成流場；transport state 拒絕後則沿用原已提供的 map |
| Abort／完整 retry | packed snapshot 相等，phase 回到 Ready、step／time 回到提交值；解除故障後可完成三個 substeps 並提交 |
| 非法重複求解 | 已完成 transport 後再求解 hydraulic／transport 會共同拒絕；接受的場、port 查詢與提交能力保留 |
| 數值參考 | flow／area／pressure／outlet 量及每個 species 分別比較 COMM_SELF，沿用 relative L2 `1e-6`／零參考 absolute L2 `1e-12` |
| 單 rank 查詢 | 只有 group rank 0 查詢 hydraulic／transport port 與 accounting，其餘 rank 等待 barrier，無隱藏 collective |
| C++ core／coupling／runtime | 全部通過；另以不含 MPI／PETSc／`-fopenmp` 的 g++ 命令編譯並執行 runtime test 通過 |
| 最終 flow graph 1+2 groups | 最大 relative L2 `7.45082e-15`，通過 |
| 最終 species graph 1+2 groups | 最大 relative L2 `3.35727e-15`，通過既有物種守恆門檻 |
| 完整 graph precommit 故障 | flow／species 各一個 2-rank 原生 CLI 作業；所有 solver rank 與 launcher 都退出 1，未 timeout、未發布 completion manifest |

precommit 案例在所有 rank 使用既有環境變數
`TUBULARFLOWIGA_INJECT_EXPLICIT_COUPLING_FAILURE_STEP=1`。這驗證完整 graph
的共同 precommit 故障與新增 adapter abort 協調相容；單 rank adapter 故障另由
上述 unit 驗證。CLI 的退出證據不等於所有 3D 內部狀態皆已量測回復。

本次未重跑原生 1D CLI／checkpoint 全套矩陣；其直接使用 runtime、沒有經過
此次變更的 adapter，前次結果見 [implicit 報告](HPC_01C_ONE_D_IMPLICIT_PROGRESS.md)。
沒有新的組裝／求解時間、加速、GPU 或跨節點證據；precommit 的 rank wall／RSS
保存在實測紀錄，不將其啟動／輸出成本當成 solver 時間。

## 證據與剩餘範圍

權威結果位於 ignored `outputs/hpc01/one-d-adapter/final/`：

- `mpi-summary.json`、1／4 threads 與 flow／species graph logs、`run_mpi.py`。
- `serial-summary.json`／log，含無 OpenMP 的完整 g++ 命令。
- `precommit-summary.json`、兩個原生 CLI 作業的各 rank 退出碼／wall／RSS／log hashes、
  `run_precommit.py` 與 `rank_worker.py`。
- 最終 build logs、`evidence-summary.json` 與 `inventory.json`，記錄 source、binary 與產物 hashes。

父目錄保留初版 19-mode 結果；它們不是最終 binary 的驗收。第一次 precommit
腳本錯將既有 `run` 匯入為不存在的 `run_rank`，solver 尚未啟動；該失敗 log、
原 worker 與 `harness_failed` summary 均保留。修正 import 後，最終版本的兩個
precommit 作業通過。所有最終 C++ 建置無新增 compiler warning／error。

HPC-01C 仍需補 3D adapters、executor 的初始配置／port routing／measurement／
before-commit／accepted-result bookkeeping、其他 runtime／CLI 邊界，以及外部
資產一致性。getter 保留本地語義不代表其呼叫者已完成 collective 邊界設計。
本協議不恢復已卡在 MPI／PETSc collective 的 rank、process loss 或 OOM-killer。
整份 [待辦清單](../WORKSTATION_HPC_TODO.md) 的其他階段保持原有範圍。
