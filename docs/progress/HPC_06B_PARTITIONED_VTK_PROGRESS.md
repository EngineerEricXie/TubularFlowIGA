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

## 接續工作

Solver 仍走既有序列輸出。下一步須從 owned elements 建立局部可視化幾何、交換所需
場值與共享點 identity；Bezier 點需沿用 extraction signature，不能以新局部編號
冒充全域共享身分。也須共同驗證 schema、cell 覆蓋及一致時間；所有 pieces 成功後才
發布 PVTU／PVD，局部失敗不得讓其他 rank 進入不匹配的 collective。

本批是格式元件，不證明大型場已免 root gather、MPI 輸出已完成或 Bezier cell 已驗收。
需以相同實際 PDE 場與既有輸出比較、測量各 rank RSS、檔案數與 metadata 成本，
再完成 HPC-06B／C。既有輸出預設與檔案介面未改。
