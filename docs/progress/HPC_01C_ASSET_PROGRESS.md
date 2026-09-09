# HPC-01C：graph 外部資產一致性

日期：2026-09-08。狀態：graph 啟動的資產檢查、原生 MPI 故障與副本數值
回歸與 CPU VCA 相容性回歸通過。HPC-01C 整體保持未勾選。

## 實作

相同設定、node／element counts 或檔案路徑，不代表各 rank 讀到相同內容。
`RunMultidomainFlow` 現在於設定快照一致後、建立 domain runtime 前，建立
邏輯資產目錄並比對檔案位元組數與 SHA-256。兩個正式 graph CLI 共用此流程。

| domain | 新增比對的輸入 |
|---|---|
| 1D | 設定引用的 SWC／OBJ network、periodic table 檔案 |
| 貼體 3D | `.ntiga`、`controlmesh.vtk`、`initial_velocityfield.txt`、periodic table 檔案 |
| 0D | 既有 model／manifest 文字快照比對維持；此次沒有新增模型格式 |

`CollectiveAssetInput.hpp` 先協調資產 key 目錄，再進行檔案讀取與逐項比對。
key 代表 domain 與輸入角色，不包含 rank 本地路徑；不同節點可使用不同本地
副本，但內容必須一致。重複邏輯 key 明確拒絕。位元組比對不會把不同的文字
格式正規化成「數值相同」，亦不取代各檔案的 parser／幾何驗證。

讀取使用固定 64 KiB buffer 與既有 SHA-256；記憶體另外包含目錄與各項摘要，
不隨最大檔案尺寸增加。POSIX descriptor 以 RAII 關閉；非阻塞 open 後檢查
regular-file 型別，避免開啟 FIFO 等候 writer。讀取期間的成長、截斷、mtime／
ctime 變化或路徑 inode 替換會拒絕。錯誤先共同協調，不讓其他 rank 越過失敗
階段開始求解。

**案例資產須在整個執行期間保持不變。** 此檢查是啟動驗證及讀取期間的變動
偵測，沒有建立檔案系統快照或鎖；後續 waveform 仍可能重新讀檔。它不能保證
偵測所有刻意修改又回復的競態，也不承諾恢復失聯的檔案系統或 MPI 程序。
本版以 Linux／WSL POSIX API 驗證；其他平台的 stat timestamp 介面尚須適配。

每個 rank 目前會額外完整掃描它可見的每份資產，包含 `.ntiga`。這是明確的
啟動 I/O 成本；大型叢集後續須量測，再評估可信 manifest／分片驗證等方案，
不能以本次小案例宣稱已具備大型共享檔案系統擴展性。

## 同時修正的本地路徑回歸

完整回歸找出先前 runtime 建構補強中的不一致：graph 入口已排除並另外驗證
應用程式路徑參數，但 runtime 再次比較 PETSc 選項時沒有沿用該清單。
因此內容相同的本地副本會在建構時被誤拒絕。

flow／transport 建構介面現在接受 caller 已驗證的 application option 排除清單；
graph 傳入原有五項：`--graph-case`、`--output-dir`、`--stop-after-step`、
`--three-d-max-newton`、`-options_file`。其他 options，包括 prefixed／unused
與已載入檔案產生的 solver 設定，仍參與比較。獨立 runtime 預設不排除任何選項。
這不改 KSP/PC 預設、物理公式、檔案格式或 trial／commit 語義。

## 驗收與限制

- `collective_asset_input_test`：world 三 ranks 與獨立一／二-rank groups，
  共 40 個故障與重試。涵蓋同大小不同內容、長度不同、key 缺少／新增／改名、
  缺檔、目錄、FIFO、device、讀取期間成長／截斷／覆寫／替換／unlink 與配置例外。
  每次確認 descriptor 數量回復，健康重試通過；三個已知 SHA 向量包含空檔、
  `abc` 與跨四個 buffer 的 196,625-byte 二進位資料。空目錄與不同本地路徑亦通過。
- `hpc_failure_regression.py`：原有 22 項加上 9 項資產案例，共 **31 項／93 份
  rank reports** 通過。有效但不同的資料庫、network、mesh、velocity、1D／3D
  waveform 各自在資產比對階段拒絕；缺檔與 FIFO 在目錄檢查階段拒絕。
  `asset-replica` 使用兩份含實際 periodic tables 的案例，正常求解並通過既有 gates。
- 正式 schema-v6 MPMD CLI：兩個 ranks 使用不同本地目錄，flow／transport 正常
  完成，五份完整 species／hydraulic 歷史 CSV 與前次驗收參考逐位元組相同。
  另以同大小不同速度檔及 rank 1 的不同 `-ksp_type` 各執行一次，均退出 1，
  沒有結果目錄。這些是小案例正確性證據，並非擴展性量測。
- 獨立 graph 群組：flow 最大 relative L2 `7.45082e-15`、species `3.35727e-15`，
  單 rank 群組為 0。既有 277 個建構故障／重試案例也通過。
- CPU VCA CLI：同版正式執行檔通過流場、傳輸、守恆與單／雙 rank checkpoint／
  restart 比較；原有數值及續跑 gate 未更動，完整輸出保留於 `vca-cli/`。

環境：GCC 11.4、OpenMPI 4.1.2、PETSc 3.15.5 real64／Int32、MUMPS，
OMP／BLAS threads 1，明確 `--map-by core --bind-to core`。單元測試外層
timeout 90 秒，CLI／群組測試 90 或 180 秒；graph harness 保留原有 child
60 秒、report rendezvous 15 秒、job 90 秒與 kill grace 5 秒。數值門檻未放寬。
CPU、OpenMP 與 graph 受影響的目標已重建，沒有新增 compiler warning。

```bash
make -C solvers/coupling collective_asset_input_test multidomain_failure_test \
  multidomain_subcommunicator_test petsc \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  timeout --kill-after=5s 90s mpiexec --map-by core --bind-to core -np 3 \
  solvers/coupling/collective_asset_input_test /tmp/asset-input-fresh
python3 scripts/hpc_failure_regression.py --output-dir /tmp/asset-graph-fresh \
  --launcher 'mpiexec --map-by core --bind-to core'
```

## 證據與接續工作

本地證據在 `outputs/hpc01/assets/`。`unit/` 保存 40 個案例，`graph-final/`
保存 31 項完整驗收及各 rank 原始退出碼／日誌；`groups/`、`construction/`、
`species-cli-replicas/`、`native-cli-negative/` 保存對應原生證據。
`source.json`／`source.tar.gz` 保存同版來源；最終彙整見 `acceptance.json`。

`graph-regression/` 保留本地路徑被建構子誤拒絕的失敗；`graph-verified/`
保留測試 fixture 漏填 periodic_table 必需 `interpolation` 欄位的失敗。
後者只修正 fixture，加入獨立 parse 檢查，未改生產 schema 或預期拒絕階段。
兩版來源也各有 archive；不把這些失敗批次列為全部通過。

後續原生 1D CLI 的資產比對亦已驗收，見
[1D 資產進度](HPC_01C_ONE_D_ASSET_PROGRESS.md)。
HPC-01C 尚待其餘獨立 CPU／explicit CLI 的外部資產檢查、adapter／
executor 局部操作與其他 runtime 邊界；浸入式仍依既有單 rank 限制處理。
本報告不宣告整份清單、HPC-01C 或跨節點 FSI 已完成。
