# HPC-06B 分片 VTK 格式元件

狀態：VTU／PVTU 格式元件與本機 ParaView 讀回通過；尚未接入 MPI solver。
日期：2026-09-10。

[PartitionedVtkOutput.hpp](../../include/PartitionedVtkOutput.hpp) 新增 `VtkPartition`、
`WriteVtuPartition` 與 `WritePvtu`。每份 piece 帶局部座標／connectivity、PointData、
CellData，以及 Int64 `GlobalPointIds`／`GlobalCellIds`。同一 piece 的 IDs 不重複，
跨 piece 允許共享點以相同 ID 出現；cell ownership 應由 collective caller 保證唯一。
不得以座標相同推論 solution-space identity 相同。

PVTU 只接受陣列 schema 與相鄰 piece 檔名，不接收場向量或完整 mesh。`GhostLevel=0`
表示不複製 ghost cells；共享邊界點仍可重複。ParaView reader 不會僅因 GlobalIds 相同
自動合併點，consumer 應保留這項區別。時間序列使用既有 PVD，每個時間指向 PVTU。
格式依據 [VTK 官方 XML 說明](https://docs.vtk.org/en/v9.6.1/vtk_file_formats/vtkxml_file_format.html)。

Writer 在開檔前檢查局部陣列長度、有限值、connectivity 範圍、完整 offsets、ID 與
schema 重複、保留名稱；PVTU 限定不重複的同目錄檔名。允許零點／零 cell 的空 piece，
但仍需帶相同空陣列 schema。正常 write／close 失敗會拋出例外。
這些檢查不是完整 VTK cell 拓撲驗證，也尚未提供 MPI agreement 或原子發布。

## 驗收

`partitioned_vtk_output_test` 產生兩個共享三角面的 tetrahedra，分成兩份非空 piece
與一份空 piece；兩個時間點 0.25／0.5。共享 GlobalPointIds 使用大於 `2^53` 的數字，
以防 reader 或日後工具經 double 轉換而遺失整數身分。另有 11 項無效輸入拒絕測試。

`test_partitioned_vtk_paraview.py` 使用本機 ParaView **5.13** 的 VTU／PVTU reader：
逐份讀回、核對空 piece、TimeValue、兩個唯一 cell IDs、cell connectivity／owner、
五個唯一 point IDs，以及全部解析速度與 scalar values；共享點座標與值完全相同。
PVD reader 的兩個時間與對應場值亦通過。讀回有八個儲存 points、五個物理 IDs，
測試不把共享點重複當成漏合併錯誤。

```bash
make -C solvers/cpu partitioned_vtk_output_test
solvers/cpu/partitioned_vtk_output_test NEW_OUTPUT
pvpython scripts/test_partitioned_vtk_paraview.py NEW_OUTPUT
```

證據：`outputs/hpc06/partitioned-v1/`，C++ exit 0，ParaView final exit 0。
最初 sandbox 中 ParaView 的 MPI singleton socket 被禁止，啟動 exit 1；同程式在有
本機 MPI 權限的環境通過，原 `paraview.log` 保留，未視為格式失敗或成功驗收。
來源、binary、fixtures 與 logs hashes 收錄於 `audit.json`。

## MPI piece 協調與索引發布

[ParallelVtkOutput.hpp](../../solvers/cpu/include/ParallelVtkOutput.hpp) 新增
`WriteParallelVtkSnapshot(comm, directory, piece, time)`：共同預檢局部資料，對
path／time／point 與 cell schema 做身分協議，root 建立不可覆寫的新 snapshot
目錄，各 rank 寫入自己的 VTU。所有 stream 正常關閉並共同確認成功後，root
才寫暫存 PVTU 並 rename 為 `snapshot.pvtu`；索引只需 O(ranks) filenames。
沒有 gather mesh／field payload。

輸出目錄必須為新路徑。失敗目錄保留 pieces 作診斷，重試使用新目錄；已發布目錄
整體拒絕覆寫。這是 live-rank 錯誤協調與索引可見性，不包含 process-loss recovery、
fsync durability、PVD 序列發布或共享檔案系統故障恢復。Global IDs／cell 唯一覆蓋
仍由 caller 建立與驗證，這層只檢查局部 IDs 與 schema。

`parallel_vtk_output_test` 在 world 3 ranks 與 split 1+2 groups 驗證：schema／time
不一致、最後 rank 的無效 tuple、真正的 rank-local 1-byte `RLIMIT_FSIZE` 寫入截斷，
以及已存在 snapshot 的拒絕。單 rank 不測不存在的跨 rank metadata mismatch。
共 13 個預期拒絕；失敗沒有 final index，換新目錄重試成功，已發布 index 與每份
piece 在覆寫拒絕前後逐位元相同。故障注入在 MPI 初始化後才設定，測後恢復。

本機 ParaView 5.13 讀回三個 groups 的 13 個成功快照，核對 point／cell ID、
各 group 獨立的座標及場值、空 rank 與 cell 數，避免混用 world communicator。
三份 rank report exit 0、無 timeout；reader exit 0。證據與 source hashes：
`outputs/hpc06/parallel-v2/audit.json`。

```bash
make -C solvers/cpu parallel_vtk_output_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 timeout --kill-after=5s 90s \
  mpiexec -np 3 solvers/cpu/parallel_vtk_output_test NEW_OUTPUT
pvpython scripts/test_parallel_vtk_paraview.py NEW_OUTPUT
```

## 共用 Bezier 點身分

`BuildBezierPointSignature` 現在抽出既有序列 builder 的 extraction 身分規則：
忽略恰為零的 coefficients，以原 1e12 quantization 與控制點 ID 建立完整排序 key，
另保留原 double coefficients 作場抽取。序列 builder 使用同一 helper，其座標一致性
檢查與共享點合併語義保留。`EncodeBezierPointSignature` 將完整 key 編為明確的
little-endian count／node／coefficient bytes，不能用單一 hash 的碰撞當成相同點。

`BezierPointOccurrenceId` 以 `64*element_id+tensor_point` 產生受 Int64 overflow
檢查的 occurrence ID；未來在相同完整 key 中選最小 occurrence，便可不依 MPI
分區取得一致代表。這是稀疏 ID，並非舊序列 mesh 的緊密 point index。
本批尚未實作跨 rank election 或代表座標／場值交換。

新 signature test 與原 `bezier_visualization_test` 都 exit 0，final build 無警告。
兩個相鄰 cubic elements 的 128 個 occurrences 得到 112 keys、16 shared points；
反轉元素與 extraction row 順序仍得到相同完整 keys／代表 IDs。相同座標但不同
solution-space ID 保持分開；另驗證 little-endian 正／負整數、Int64 邊界及八項
錯誤拒絕。重現：`make -C solvers/cpu bezier_point_signature_test bezier_visualization_test`，
再執行兩個 executable。來源與 logs 見 `outputs/hpc06/signature-v1/audit.json`。

## 分散式完整 key 代表選定

[DistributedPointIdentity.hpp](../../solvers/cpu/include/DistributedPointIdentity.hpp) 新增
`ResolveDistributedPointIdentities`。先按 occurrence ID 分配檢查全域唯一性，再按 key
hash 分配完整 key，以完整 bytes 比較、選最小 occurrence，回傳該 ID 及其來源 rank。
Hash 只決定負責配對的 rank，沒有用 hash 相同取代完整 key equality。來源 rank 可隨
分區改變，代表 occurrence ID 保持相同；回傳順序與 caller 的局部輸入一致。

每個階段的局部配置／解析錯誤共同協調，使用 Alltoall／Alltoallv 而非全目錄 gather。
每 rank 各次交換的送出／接收 wire payload 均設上限，預設 64 MiB，局部 occurrences
上限一百萬，並限制 MPI int counts。分配偏斜可觸發接收端 cap，全群拒絕後仍可重試。
這些是 wire／輸入數量限制，不是 RSS 上限；buckets、flat buffers 與 key map 的
暫存可能同時存在。沒有 process-loss recovery 或跨節點量測宣稱。

實際 MPI 測試在 3／1／2-rank groups 各跑三種分區（全放 rank 0、交錯、反向交錯），
並反轉局部輸入順序。131 occurrences 包含兩個 cubic 元素的 128 點與三個 binary-key
案例，得到 114 個完整 keys；核對每筆最小 ID 及正確 communicator owner。
另含全空群組、內嵌 NUL、相同前綴不同完整 key、大於 2^53 的 Int64 ID。
14 項預期拒絕涵蓋重複 occurrence、空 key、局部數量配置、送出 cap 與接收偏斜 cap，
每項核對指定錯誤階段並成功重試。三份 rank reports exit 0、無 timeout。

```bash
make -C solvers/cpu distributed_point_identity_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 timeout --kill-after=5s 90s \
  mpiexec -np 3 solvers/cpu/distributed_point_identity_test
```

證據與 source／binary／logs hashes：`outputs/hpc06/identity-v2/audit.json`。
此元件只交換身分 metadata；代表點座標、實際場值的交換與局部 Bezier piece 建立
尚未接入，故不能宣稱 solver 分片輸出已完成。

## 代表點資料交換

[DistributedPointValues.hpp](../../solvers/cpu/include/DistributedPointValues.hpp) 新增
`ExchangePointRepresentativeValues`：以代表 occurrence／rank 發送局部請求，owner
查找其局部 row，回傳共同 component 數目的 tuple。座標與多個場可放在同一 tuple。
資料以 double bits 的 little-endian 編碼交換，回傳順序與局部 occurrence 對應；
本機代表也經相同請求檢查。來源 values 不修改，未 gather 完整場到 root。

共同預檢 components／形狀、owner、局部 ID、finite values；對缺失代表、錯誤回應
row／來源與 wire cap 共同拒絕。沿用身分元件的每次送／收 wire caps，不宣稱為
allocator peak 上限。這層相信 caller 提供已驗證的代表配對；不重新比對 extraction key。

同一 `distributed_point_identity_test` 在 3／1／2-rank 三種分區中檢查每筆交換值的
代表 ID、owner rank 與 signed zero，並測試全空群組。新增 16 項指定階段拒絕：
無效 owner、缺失代表、非有限來源、跨 rank component 不一致、接收偏斜 cap，
以及 owner 的回應建構 cap。每項核對來源 values 位元不變並成功重試；原 14 項
identity 拒絕與全部 ID 一致性仍通過。三份 rank reports exit 0、無 timeout。
命令同上一節，證據為 `outputs/hpc06/point-values-v2/audit.json`。

這是 tuple 交換元件的驗收，還須用 solver 的局部 extraction 建立真正的座標與場
資料、檢查共享座標一致性，再組裝 Bezier piece 與發佈時間序列。沒有宣稱已完成
實際 PDE 分片輸出或大型 RSS 驗收。

## Owned-element Bezier 分片建構與內插修正

[ParallelBezierVisualization.hpp](../../solvers/cpu/include/ParallelBezierVisualization.hpp) 新增
`BuildParallelBezierPartition`：caller 提供該 rank 的 owned elements、全域元素數、
field schema 與只存取本地 owned／ghost 控制值的 callback。每個 occurrence 建立
完整 signature、座標與 extracted field tuple，選代表並交換後，僅合併該 piece
內相同代表 ID。輸出有 GlobalPointIds、GlobalCellIds、owner 與 HigherOrderDegrees。

共同驗證 owner／元素範圍、全域元素總數、occurrence 唯一性；三者合併證明全元素
恰好出現一次。全域僅交換每 rank cell count 及 bounds 等 metadata。共享座標沿用
序列 tolerance，canonical 座標合併後，每 cell 的 64 Gauss locations 必須有正且
有限的 Jacobian，並拒絕 collapsed point IDs。此元件沒有重做跨元素 volume overlap
認證；該認證仍須由上游幾何驗證提供，不能以局部 Jacobian 替代。

在 world 3 ranks 與 split 1+2 groups，各三種 owned-element 分區並反轉本地元素順序，
九個快照皆與原序列 builder 的座標／速度／scalar 逐位元相同。15 項指定拒絕涵蓋
缺元素、重複元素、共享座標衝突、缺控制值與負 Jacobian。所有 rank exit 0。
ParaView 5.13 讀回九個真正的 cubic Bezier 快照，檢查 112 個共享點 IDs、兩個 cell
IDs、degree (3,3,3)、完整場，以及 cell 內參數點 (0.2,0.4,0.6) 的座標／場內插。
解析速度對座標的最大差為 2.220446049250313e-16；序列／分片比較仍是逐位元相同。

這個 reader 測試揭露兩項排序問題，已修正：

- 共用 `VtkCubicHexTensorIndices` 原本轉置兩個面的內部點順序。以本機 VTK
  `PointIndexFromIJK` 核對後修正四個 entries。重建舊排序 VTKHDF 的 unit cube
  在該參數點實際讀出 (0.2,0.4179712,0.5820288)，最大座標誤差約 0.0179712。
- 新分片 XML 若宣告 0.1，reader 會採舊版高階 connectivity 相容轉換，額外交換兩條
  邊的點。以同一 v2 fixture 的獨立副本只改版本為 2.2，九個內插測試即通過。
  正式 partitioned VTU／PVTU writer 現在輸出 XML 2.2，v3 的原生檔案亦全部通過。

既有低階 `WriteVtu` 介面不變。原 `bezier_visualization_test`、VTKHDF schema／resume
與基本分片 VTU/PVTU reader 均通過。`vtkhdf-paraview-test` 新增 unit-cube 的元素內插
驗證；修正後通過。舊排序的 VTKHDF 因 GeometryHash 不同而拒絕續寫，測試已確認；
需要重新輸出到新檔案，不能把新時間幀接在舊 connectivity 後。數值 checkpoint 格式
與求解器積分不受此可視化排序修正影響。

失敗證據保留於 v1／v2：最初解析測試錯要求 velocity 與 coordinate 位元相等
（1+2/3 與 5/3 有 rounding 差），改為明確 1e-14 解析門檻後才進一步抓到真正的
內插誤差。舊 HDF 的首次 reader 驗證停在單幀檔案無 TimestepValues，因此另用
vtkHDFReader 直接量測上述內插誤差；未把前者算作幾何失敗證據。

最終證據：`outputs/hpc06/parallel-bezier-v3/audit.json`，前輪分別保留於
`parallel-bezier-v1`／`parallel-bezier-v2`。重現建構與 reader：

```bash
make -C solvers/cpu parallel_bezier_visualization_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 timeout --kill-after=5s 120s \
  mpiexec -np 3 solvers/cpu/parallel_bezier_visualization_test NEW_OUTPUT
pvpython scripts/test_parallel_bezier_paraview.py NEW_OUTPUT
```

尚須由真正 PETSc owned／ghost state 提供 callback，接入 CLI／runtime 輸出頻率與
PVD 時間序列，再驗收實際 PDE 數值及大型場 RSS；本批不是完整 solver 整合。

## 接續工作

Solver 仍走既有序列輸出。局部 Bezier 建構、shared IDs／代表 tuple 交換、cell 覆蓋
與 PVTU 協調已有上述元件；接下來須接入 PETSc owned／ghost state、CLI 與 PVD
時間序列。完整幾何 overlap 認證的來源與生命週期也須明確，不能以局部正 Jacobian
代替跨元素檢查。

需以相同實際 PDE 場與既有輸出比較，測量各 rank RSS、檔案數與 metadata 成本，
再完成 HPC-06B／C。尚不證明大型 solver 場已免 root gather。既有預設選項保持，
但上述可視化排序修正會拒絕沿用舊排序 VTKHDF 的 geometry hash。
