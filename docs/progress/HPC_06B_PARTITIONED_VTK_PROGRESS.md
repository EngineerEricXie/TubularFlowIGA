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

## 接續工作

Solver 仍走既有序列輸出。下一步須從 owned elements 建立局部可視化幾何、交換所需
場值與共享點 identity；Bezier 點需沿用 extraction signature，不能以新局部編號
冒充全域共享身分。schema／時間協議與 PVTU 發布已有上述元件；尚須整合 cell 覆蓋驗證、PVD
序列發布及 solver 的局部場交換，不可只因格式層通過就宣稱整個流程完成。

本批是格式元件，不證明大型場已免 root gather、MPI 輸出已完成或 Bezier cell 已驗收。
需以相同實際 PDE 場與既有輸出比較、測量各 rank RSS、檔案數與 metadata 成本，
再完成 HPC-06B／C。既有輸出預設與檔案介面未改。
