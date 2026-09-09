# HPC-01C：配置一致性與建構前檢查

狀態：**本次列出的驗收通過；HPC-01C 整體仍部分完成。**
接續 [首批一致失敗處理](HPC_01C_PROGRESS.md)，處理 rank 各自讀到合法輸入、
但設定不同，導致後續 collective 順序或矩陣尺寸不一致的路徑。

## 已實作的保證

1. `RequireCollectiveSameText` 使用 communicator rank 0 的內容逐段比較。
   64-bit 長度與 4096-byte 固定 buffer 避免為廣播再複製整份 root 文件。
   即使長度不符、內容不同或文字含 NUL，所有 rank 都完成同一串廣播後，
   才透過共同失敗協議退出。測試包含跨兩個 chunk 後才出現的差異。
2. graph 在集體 runtime 建構前，比較解析後的 stop-after-step、Newton 上限
   與 failure-injection step，以及實際供 parser 使用的 graph manifest、
   0D model、native 1D／貼體 3D JSON 快照。超出 graph horizon 的停止步數
   也提前在 collective local stage 拒絕。
3. `ReadMultidomainConfiguration` 與 `ReadZeroDFlowModelConfiguration` 增加
   回傳輸入快照的 overload；原有單參數 overload 保留。0D 模型的快照是
   實際解析的字串，不是解析後重新開檔做一次不相干的比較。1D／3D runtime
   的前置 parser 也直接使用已比較的快照。
4. 只比較內容與 logical domain identity，不把絕對路徑作為相等條件。
   相同的相對路徑配置與資產可複製到不同本地目錄，三 rank replica case
   已通過原生流量／守恆 gate。兩個獨立 communicator 可以使用不同 graph。
5. `OwnedRowAssembler` 在建立 PETSc 物件前，額外確認全域 node、element、
   field counts 一致。正數但互異的 field count 也會全員拒絕，不只拒絕零／負值。
6. transport 的配置／label 複製、coupling pattern／element scratch 建立、
   volume-rule 綁定、boundary 與 initial-value 解析，移至共同 local stage，
   並在 PETSc 物件建立前完成。原本在 member initializer 丟出的錯誤 field
   index 現在能一致傳回。flow 也先協調完整 boundary array shape、constraint
   mask 複製及 owned-element 載入，再進入 boundary catalog 通訊。

JSON 比較是**位元組一致性**，不是語義正規化。各 rank 的空白、鍵順序或
字串表示不同，也會拒絕；單 rank 的 JSON 接受規則不變。此規則讓啟動時的
輸入副本有可檢查的身分，沒有聲稱每一種等價 JSON 都可混用。

本次不修改 `.ntiga`、cache、partition、velocity、checkpoint 格式，不改
元素公式、時間積分方法、SI ports、trial／commit／rollback 或數值 gate。

## 驗收與環境

沿用 WSL 工作站、GCC 11.4.0、Open MPI 4.1.2、PETSc 3.15.5 real/double、
32-bit PetscInt。MPI 測試使用 3 ranks，每 rank OMP／OpenBLAS threads 為 1。
基準 HEAD 為 `ee8da2ab528849814b6d12c8186476b5f7f59124`，但工作樹包含未提交修改；
精確來源、執行檔與輸入身分以本次 `inventory.json` 為準。

```bash
make -C solvers/cpu collective_failure_test vca_3d_runtime_test \
  parallel_ownership_test body_fitted_subcommunicator_test petsc \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
make -C solvers/coupling multidomain_failure_test multidomain_subcommunicator_test \
  multidomain_config_test petsc \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
python3 scripts/hpc_failure_regression.py --output-dir /tmp/hpc-agreement-new-run
```

輸出目錄必須尚不存在。harness 沿用 child 60 秒、report rendezvous 15 秒、
job 90 秒與 kill grace 5 秒的限制。它清除繼承的兩個 coupling failure-step
環境變數，再由指定案例只在 rank 1 設定；正常案例不受外部測試 knob 影響。
MUMPS direct options 用於 harness 與 subcommunicator 參考比較；VCA 原生
單 rank 測試仍使用其預設迭代求解器。

`outputs/hpc01/agreement/verified/summary.json` 的 **15 個案例全部通過**：

| 案例 | 預期／實際 MPI 作業退出碼 | 驗證 |
|---|---:|---|
| unit | 0 | 共用錯誤協議、text agreement、零 owned-row rank、建構／組裝拒絕、rollback |
| arguments、input、database、output | 1 | 首批四種 production 故障，指定階段診斷與無成功標記 |
| stop、newton、injection | 1 | 只有 rank 1 的合法控制值不同，建構 collective runtime 前退出 |
| manifest | 1 | rank 1 使用合法 Aitken graph，其餘使用 explicit graph |
| one-d-config | 1 | rank 1 的 cells-per-segment 不同，兩份配置均可獨立解析 |
| three-d-config | 1 | rank 1 的 3D dt 不同，兩份配置均可獨立解析 |
| zero-d-config | 1 | rank 1 的 source 初始壓力不同，模型與 graph 均可解析 |
| healthy、replica | 0 | 非零 flow graph 與不同路徑的相同輸入副本，原生物理 gate |
| zero-d-healthy | 0 | 8-step 0D／3D clock、source／RCR history、accounting 與 manifest gate |

所有 fault case 均為正常 PETSc finalize 後回傳 1；不是 timeout 或以 MPI_Abort
代替受控失敗。wrapper 的原始非零 `run.json` 保持 failed，另由測試 summary
判定是否符合故障契約。每項均驗證三份 rank reports、log hash 與確認訊息。

這批微型 case 最大 rank wall 為 `0.714662 s`（0D healthy），最大單 rank
RSS 為 `45,395,968 bytes`（replica）。包含啟動與量測成本，不含 report rendezvous；
不是 isolated performance repetitions，也不是 scaling 證據。

## 數值與相容性回歸

- 原有 `multidomain_config_test` 通過，保留 schema／0D 角色接受與拒絕規則。
- `vca_3d_runtime_test` 通過，涵蓋 native flow／transport、boundary、trial、
  commit／abort、outlet 模型與狀態回復。
- `multidomain_subcommunicator_test` 的 1+2 程序群均通過；flow graph 對串行
  最大相對 L2 為 `7.45082e-15`，species graph 為 `3.35727e-15`，
  各單 rank 群組都是 0。數值、edge identity、迭代次數與物種 accounting gate 保留。
- 1005-node／720-element 貼體案例的 `body_fitted_subcommunicator_test` 通過。
  兩 rank 群組相對誤差：velocity `3.48315e-15`、pressure `9.94533e-16`、
  transport `1.70503e-14`、mass `1.62206e-14`；單 rank 群組都是 0。
  同時驗證正 Jacobian、halo 更新、獨立 ownership catalogs、trial rollback／retry
  與 checkpoint readback。三 rank 均退出 0，無 timeout；最大 wall `81.243808 s`，
  最大單 rank peak RSS `266,563,584 bytes`。原有 CPU relative L2 `1e-6` 不變。
- 受影響 CPU／coupling CLI 與測試完成重建，無新增 C++ compiler warning。
  未重跑整套耗時較長的 multidomain CLI suite；不將舊完整 suite 紀錄套用到此版。

本次是正確性與失敗處理工作，沒有宣稱速度改善。組裝／solve 的既有 runtime
log 保留；分階段效能歸因、獨立重複統計、CUDA 與多節點驗收不屬於這次證據。
固定 buffer 比較仍增加啟動通訊，後續 HPC-00C／09D 必須量測整體成本。

## 證據與仍待完成的範圍

本地產物位於 ignored 目錄 `outputs/hpc01/agreement/`：
`verified/summary.json`、`subcommunicator.json`、`body/ranks/rank-*/run.json`、
`regression-summary.json`、建置／回歸 logs 與 `inventory.json`。
原始 `initial` 執行保留為建構補強前的結果；以 `verified` 與本次 binary hashes
辨識最後版本。首批測試的缺 report／launcher status 問題仍保存在前一份報告的目錄。

HPC-01C 尚不可勾選，剩餘工作包括：

1. graph 啟動的可見 PETSc 選項已在後續工作中加入比較，見
   [PETSc 選項驗收](HPC_01C_PETSC_OPTIONS_PROGRESS.md)。本文的 15-case
   結果保留為先前版本證據；其他 CLI 的選項檢查仍待補齊，啟動快照也不代表
   後續選項修改或預先建立的 PETSc 物件具有相同設定。
2. 相同配置內引用的 SWC／geometry／waveform／`.ntiga` 內容仍可能不同。
   node／element／field counts 只是 shape，沒有取代內容、field-name/order、
   partition ownership 與物理身分檢查；不以此宣稱完整 asset identity 驗收。
   後續 graph 啟動已加入外部資產內容比對，見
   [資產一致性驗收](HPC_01C_ASSET_PROGRESS.md)；其餘獨立 CLI 仍待補齊。
3. caller-side 參數配置、其他 constructor initializer、halo／boundary catalog
   配置、PETSc 部分建立後清理，以及其他 CPU／1D／explicit CLI 的局部錯誤邊界。
4. 後續 port 更新、coupling bookkeeping、診斷、viewer／checkpoint／輸出與
   PETSc 回傳碼處理；新增 runtime 時亦須接入協議。
5. `CollectiveLocalStage` 仍要求同一 communicator 的全員以同樣順序呼叫，
   callback 不得含 collective。`RequireCollectiveSameText` 本身是 collective，
   不可放進其 local callback。此工作不提供 process-loss、OOM killer 或節點失聯恢復。

原始完整 goal 仍依 [待辦清單](../WORKSTATION_HPC_TODO.md) 持續追蹤；不縮減為
上述已通過的微型案例或已完成的子任務。
