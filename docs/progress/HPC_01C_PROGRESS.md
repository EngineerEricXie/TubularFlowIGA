# HPC-01C：一致失敗處理進度

狀態：**部分完成，HPC-01C 保持未勾選。** 本報告只涵蓋下列已接入的
失敗邊界，不宣稱所有 runtime、CLI 或 PETSc 內部錯誤均可一致退出。
本文保留首批實作與當時六項測試的證據；後續配置一致性與建構階段補強，見
[配置／建構補強進度](HPC_01C_CONFIGURATION_PROGRESS.md)。

## 本次實作

- [CollectiveFailure.hpp](../../solvers/cpu/include/CollectiveFailure.hpp) 的
  `CollectiveLocalStage` 在本地 callback 完成或捕捉例外後，以 supplied
  communicator 選出最小失敗 rank，再廣播相同的階段、rank 與錯誤訊息。
  捕捉 `std::exception` 與未知例外；1024-byte stack buffer 避免在錯誤協議前
  再配置診斷字串，過長訊息截斷。協議通訊納入既有 Communication 計時。
- `RunMultidomainFlow` 的參數解析、graph 輸入、domain preflight 及 root
  輸出使用共同協議，包含原本在 preflight 外的參數解析。
- `OwnedRowAssembler` 建構時的 required-element 讀取／範圍檢查，以及
  dense field pattern、矩陣 adjacency／預配置陣列建立，先完成本地錯誤協調
  才進入 PETSc 矩陣建立。**建構與 CreateMatrix 要由 communicator 全員呼叫。**
- 貼體 flow 的元素積分／插入、continuity 與 boundary-flow 本地積分、
  boundary values／插入，以及 transport 的 trial 輸入、boundary 解析、
  元素積分／插入與 boundary 插入，均在後續 collective 前協調本地例外。
- [PetscReadArray.hpp](../../solvers/cpu/include/PetscReadArray.hpp) 在例外展開時
  還原 flow 的 current／history／owned／residual read view，包含零長度 Vec。
  已失敗的 trial 仍不可 PrepareCommit，保留原有 AbortStep／rollback 語義。

不改元素公式、邊界定義、數值容許值、檔案格式或 communicator 擁有方式。
新增同步有通訊成本；這批測試未量測 scaling 或宣稱效能改善。

## 協議前提與尚待完成項目

callback **只能執行本地操作**。所有 rank 必須呼叫同樣階段且順序相同。
將一個含 collective 的大函式整體包在 callback，不能解決其他 rank 已進入
該 collective 的等待。此 helper 也不是 communicator 階段／配置一致性驗證器。

尚待補齊：

1. graph manifest、0D／1D／3D JSON、停止步數、Newton 上限及 failure knob
   已補上快照一致性檢查。graph 啟動的可見 PETSc 選項亦已
   [完成獨立驗收](HPC_01C_PETSC_OPTIONS_PROGRESS.md)。仍需檢查外部 geometry／
   waveform／packed database 身分與其他 embedding／CLI 配置入口。
2. 其餘 CPU／1D／explicit CLI 的資料讀取與局部組裝邊界；graph 後續 runtime
   建構、halo／boundary catalog 配置、port 更新、耦合 bookkeeping 與診斷。
3. transport 的大型配置複製／coupling pattern 建立已移入共用 local stage，
   flow 初始 boundary shape／owned-element 載入亦已協調。仍需補齊 caller-side
   參數配置、其他 member initializer，以及 PETSc 資源部分建立後的清理契約。
   只捕捉 `bad_alloc` 的測試不代表實際記憶體耗盡可完整恢復。
4. 補齊 PETSc 回傳錯誤檢查與既有 viewer／checkpoint／輸出路徑；本報告的
   graph output gate 不代表所有分散式 I/O 路徑均完成。
5. 更完整的故障矩陣及整合驗收，才可勾選 HPC-01C。

本協議依賴每個 MPI 程序仍存活並可進入通訊；不承諾 SIGKILL、OOM killer、
節點失聯、MPI 內部故障後原地續跑，也沒有以 MPI_Abort 通過受控故障案例。
測試程式只有遇到**不符合預期的測試失敗**才用 MPI_Abort 結束，以保留非零證據。

## 可重現命令

環境沿用 HPC-01A/B：WSL 工作站、GCC 11.4.0、Open MPI 4.1.2、
PETSc 3.15.5 real/double、32-bit PetscInt，單 rank 一個 OpenMP／OpenBLAS thread。
基準 HEAD 為 `ee8da2ab528849814b6d12c8186476b5f7f59124`；本次含未提交修改，
應以此次 inventory 的 working-tree／binary hashes 識別，不能僅以 HEAD 當作來源身分。

```bash
make -C solvers/cpu collective_failure_test vca_3d_runtime_test parallel_ownership_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
make -C solvers/coupling multidomain_failure_test multidomain_subcommunicator_test petsc \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
python3 scripts/hpc_failure_regression.py --output-dir /tmp/hpc-failure-new-run
```

輸出目錄必須不存在。harness 預設 Open MPI launcher `mpiexec --oversubscribe`，
可用 `--launcher` 指定其他安裝；它在每 rank 使用既有 `hpc_rank_run.py` 記錄
stdout／stderr、hash、退出碼及 GNU time RSS。每 rank child timeout 60 秒，
measurement wrapper 另以最多 15 秒等待其他 rank 的 report 完成，再回傳原本
退出碼，避免 Open MPI 因第一個非零退出提早殺掉尚在寫 report 的 wrapper。
此測試收集協議需要所有 rank 可見同一個新建輸出目錄，不是 runtime 的故障復原機制。
作業 timeout 90 秒、額外 kill grace 5 秒。故障作業必須真的退出 1，不能是 0、
124、signal 或缺少 rank 紀錄；原始 `run.json` 的失敗狀態不會改寫成成功。

MPI socket 在 sandbox 中受限；數值與 MPI 驗收在允許本機 MPI 通訊的執行環境進行。
上述微型案例不是 cluster 擴展性測試；大型案例仍需排程分配的計算資源。

## 驗收範圍

| 案例 | 故障與驗證 |
|---|---|
| 共用協議 | 三 rank、1+2 子 communicator；標準例外、bad_alloc、未知例外、多 rank 同時失敗與過長訊息；健康階段可繼續 |
| Owned rows | 只有 rank 1 的 field count／matrix pattern 錯誤；rank 2 為零 owned rows 的 partition owner |
| 傳輸輸入 | 只有 rank 1 的 velocity catalog 長度錯誤；AbortStep 精確恢復 committed state，重試符合原先 CPU relative L2 ≤ 1e-6 |
| 傳輸組裝 | 只有 rank 1 的 compiled dt 為零，實際元素 kernel 丟錯；所有 rank 不得 commit 且 AbortStep 還原 |
| 流場組裝 | 只有 rank 1 的 density 非法，實際 Navier–Stokes kernel 丟錯；拒絕 commit、還原狀態後重做同一失敗 trial，驗證 read views 已釋放 |
| Graph arguments | 只有 rank 1 帶未識別參數；全 rank 回傳 1，rank 0 顯示 rank 1 診斷 |
| Graph input | 只有 rank 1 指向不存在的 graph root；全 rank 回傳 1 |
| Graph domain preflight | rank 1 使用存在但截斷的 `.ntiga`；確認錯誤來自 domain preflight |
| Graph output | 輸出目錄的 parent 是 regular file；真的執行求解後在 root writer 失敗，全 rank 回傳 1 |
| Healthy graph | 三 rank 執行既有四 domain fixture，完成 native flow/conservation 檢查及成功標記 |

受控 graph 失敗均確認沒有 `graph_binding_manifest.json`。harness 另外檢查
失敗訊息的階段，避免將「測試沒進入預定路徑」誤認為成功故障注入。
小型 algebraic fixture 用於所有權、同步與狀態恢復，不當作完整物理驗收。
非零 flow／species graph 的數值相容性由原生 fixture 與串行比較另外驗證。

## 首批結果與證據（當時版本）

首批六項結果位於 `outputs/hpc01/failure/verified-reports/summary.json`，
保留被測 `collective_failure_test`／`multidomain_failure_test` 的 SHA-256。
每項都驗證三份 rank report、log hashes、未逾時、退出碼與 case-specific
確認訊息；四個故障案例的 launcher／全部 child 退出碼都是 1。

| 案例 | 作業退出碼 | 最大 rank wall（秒） | 最大 rank peak RSS（bytes） |
|---|---:|---:|---:|
| unit／runtime 局部故障 | 0 | 0.314070 | 40,591,360 |
| arguments | 1 | 0.314009 | 33,275,904 |
| input | 1 | 0.313992 | 34,689,024 |
| database | 1 | 0.314043 | 35,500,032 |
| output | 1 | 0.414466 | 44,126,208 |
| healthy | 0 | 0.414225 | 44,052,480 |

wall 包含 wrapper 測量／process startup，不包含後續 report rendezvous；
不是 assembly 或 solver time。這不是效能工作或獨立重複量測，不從這些數字推論
加速；本次分階段效能、CUDA peak、跨節點量測均 N/A。

另已通過：

- `multidomain_subcommunicator_test`：1+2 群組的 flow graph，對各自串行解的
  maximum relative L2 為 `7.45082e-15`；species graph 為 `3.35727e-15`。
  group 0 都是 0。沿用 CPU `1e-6`、解析零量 absolute `1e-12` gate，
  包含原生守恆／edge／domain／global species accounting 檢查。
  命令、log 與 binary hashes 在 `outputs/hpc01/failure/subcommunicator-verified.json`。
- `vca_3d_runtime_test`：既有單 rank flow／transport／boundary／trial 回歸通過，
  此測試沿用預設迭代求解器。最終 log 為 `vca-verified.log`。
- `parallel_ownership_test`：三 rank 的共享元素、獨立 partition／row catalogs、
  empty rank 與精確覆蓋故障測試通過，log 為 `failure-ownership.log`。
- CPU `petsc`、coupling `petsc` 及上述新增／受影響測試均完成重建，無新增 C++
  compiler warning。編譯 wrapper 的 sandbox socket 提示不當作 MPI runtime 驗證。
- 未重跑耗時較長的完整 multidomain CLI 全模式 suite；HPC-01A 的舊完整紀錄
  不能當成此版的完整回歸。此次以實際 runner entry point 的故障、三 rank
  正常 graph、flow／species 子 communicator 及 VCA 原生測試作為明確限定的證據。

上述結果與建置 logs 保存在 `outputs/hpc01/failure/`；此目錄為 ignored
本地測試產物，不提交二進位或生成案例。`inventory.json` 記錄本次工作樹、
執行檔及證據身分。完整 hardware／PETSc 設定仍參照 HPC-00A／01A/B inventories。

## 已知測試修正紀錄

- 初版 transport retry 使用逐位元比較，新 GMRES 求解只有 `2.22045e-16`
  相對 L2 差異。改用預先登記的 CPU `1e-6`；沒有放寬既有數值門檻。
  rollback 的精確比較仍保留。原始失敗 log 保留於本次證據目錄。
- harness 初版設定 Open MPI `orte_abort_on_non_zero_status=0`，實測讓 launcher
  回傳 0，即使各 rank 回傳 1。驗證器因此拒絕該次執行；已移除此設定並以
  預設非零退出傳遞重新驗證。原始 `attempt-01` 保留，不作為通過證據。
- 預設非零退出下曾遇到 rank 1 的 child 已結束，但 measurement wrapper 在
  寫入 `run.json` 前被 launcher 結束；該次 `verified` 因缺少證據被拒絕。
  已加入上面的 bounded report rendezvous，仍保留真實作業非零退出，沒有將
  測試 child 的失敗改寫成 0。較早通過的 `attempt-02`／`final` 也保留，
  最終來源與結果以本次報告指定的最後一組 binary hashes 為準。

完整清單仍由 [WORKSTATION_HPC_TODO.md](../WORKSTATION_HPC_TODO.md) 追蹤。

後續的貼體 3D port 查詢、flow adapter 輸入／求解準備及相關 MPI 故障驗收，
見 [3D port 進度](HPC_01C_THREE_D_PORT_PROGRESS.md)；該報告另列尚未補齊的
3D runtime／staged adapter 邊界，不改變 HPC-01C 仍部分完成的狀態。

其後的 runtime trial／Abort／快照回復與 outlet evaluation 補強，見
[3D trial 進度](HPC_01C_THREE_D_TRIAL_PROGRESS.md)，包含新的 48 個故障案例
與既有 port／graph／VCA 回歸。

staged 輸入發布、傳輸邊界、物種收支與 mass／source reduction 的後續補強，見
[3D staged 進度](HPC_01C_THREE_D_STAGED_PROGRESS.md)，72 個故障案例及流向反轉驗收通過。

其後的 ReferenceBoundaryFlow／Summary 共同查詢補強、64 個 port／診斷案例
與同版回歸，見 [流場診斷進度](HPC_01C_FLOW_DIAGNOSTICS_PROGRESS.md)。

Newton assembly／scatter／KSP setup／solve／update、收斂決策與 logging 的
後續協調，見 [Newton 進度](HPC_01C_NEWTON_PROGRESS.md)，75 個 trial 案例、
同版 graph／VCA 回歸及 profile 時間分解檢查通過。

flow 初始化候選發布、resolved 邊界比較與 transport 全域 gather 的補強，見
[初始化進度](HPC_01C_INITIALIZATION_PROGRESS.md)，93 個 trial、83 個 staged、
12 個 gather 故障案例及同版回歸通過。
