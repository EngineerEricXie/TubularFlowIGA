# HPC-07 分散式 FSI 接線契約與缺口

狀態：原始碼契約盤點完成，HPC-07A–D 尚未驗收。
基準：`17d620c`，2026-09-10。這份盤點不啟動求解，以避免干擾正在執行的
固定網格 rank 量測。數值／效能／CUDA allocation 對本次原始碼盤點為 N/A。

## 現有能力與實際限制

| 元件 | 已有契約 | 尚缺的分散式行為 |
|---|---|---|
| `DistributedSurfaceInterface.hpp` | 不可變 material node ID、owned 發布權、reference／layout／partition hashes、reference lumped areas、kinematics／traction stamp | reference positions／triangles 完整複製；拒絕 owned node 為空；沒有 ghost 通訊與三角形貢獻所有權 |
| `SurfaceOwnershipValidation.hpp` | communicator rank／layout 一致、root 提供 reference ID catalog，sharded incidence 驗證節點唯一發布與覆蓋 | 不分散 reference 幾何、surface quadrature、力或膜求解；不能當成完整 FSI 驗收 |
| `MaterialSurfacePatchMap.hpp` | reference material／patch／triangle provenance 綁定 | 建構要求完整單分區 |
| `FluidSurfaceTraction.hpp` | 原 Cauchy traction、barycentric consistent force／projection、合力與力矩驗證 | 要求完整單分區；遍歷所需狀態與 bounded consistent projection，未路由跨 partition 的 nodal contributions |
| `MovingImmersedTransientFlowFsiRuntime.hpp` | moving epoch、traction／kinematics stamp 與 trial commit | fluid／structure layout 均須單分區；capture 走完整 geometry／state 假設 |
| `PretensionedMembraneFsiRuntime.hpp` | 膜 model、reference normals、clamped state 與 trial rollback | wrapper 要求完整單分區 |
| `StrongFluidStructureCoupling.hpp` | Aitken、收斂、prepare／commit／abort、結果 hashes | 明確限定 partition_count=1；AreaTotal 是 local sum，不能直接作 global residual denominator |
| `DynamicWeightedAitkenRelaxation.hpp` | 接受 externally reduced numerator／denominator／scale，分離 globally comparable control identity | 尚未在正式 FSI 接上 MPI reduction；local weight total 必須正，拒絕空 rank |

## 實作時必須維持的分區契約

1. **Material node 發布權。** 全域每個 reference node 恰一個 owner；ghost 僅為
   本地 geometry／quadrature／structure stencil 所需資料，不重複發布或計入
   residual。rank 可以沒有 owned nodes 或沒有 owned fluid surface quadrature，
   仍參與相同 collective sequence。先保留既有 replicated reference 相容路徑，
   明列其 O(N+T) 每 rank 成本；新增局部表示時不得假裝原 layout 已免複製。
2. **Reference 面積。** Immutable area 權重來自 reference triangles，各 triangle
   恰一次貢獻 A/3 到三個 material node，路由並加總到 node owner。三角形跨越
   node 分區並不使它在每個 node owner 各積分一次。需驗證全域 nodal area sum
   等於 reference triangle area sum，ghost 不計入總和。
3. **流體牽引力。** 積分 authority 由 owned fluid cell／surface quadrature 決定，
   不由 nodal owner 決定。同一 material triangle 可穿越多個 fluid cells；其
   clipped quadrature contributions 不能因 triangle ID 相同而誤刪。沿用原
   canonical triangle／barycentric provenance 與 `FluidOnStructureCauchyTraction`
   的方向約定。每個積分點的 Ni·traction·weight 路由到 material node owner；
   consistent projection 的 mass matrix contributions 也須使用同一 authority。
4. **場交換。** Owner 接收加總後的 consistent nodal force／traction，displacement
   與 velocity 由指定 structure authority 發布並送到需要的 ghost。Time、step、
   coupling iteration、reference／layout／partition／producer state stamps 都需
   檢查；不能將全域 topology hash 當作 partition identity。
5. **集中小膜的範圍。** 第一版可由單 owner 解小膜及必要的 bounded projection，
   但 fluid traction quadrature 必須真的分散；記錄集中配置大小、通訊、RSS 與
   solve time。不能把每 rank 完整重算膜稱為 distributed structure solve。
6. **強耦合。** 面積加權 RMS numerator／denominator 使用 owned nodes 的全域和；
   最大殘差與 scale 用全域 max。Aitken numerator／denominator 亦需全域 reduction，
   空 rank 貢獻 0。所有 ranks 使用同一 proposal／接受決策。單 rank 例外必須先
   共同傳播，再一起 rollback／retry，不能有些 ranks commit、有些 ranks abort。
7. **續跑。** 接受的 structure state、moving geometry、donor／Aitken 契約依
   HPC-05 保存或明確重建。不同資源映射的處理不能只改 partition stamp 來繞過
   身分檢查；需另訂合法重新分區及全域 material ID 對應。

## 必要驗收與實作順序

先做 HPC-07A 的空分區／node ownership／reference triangle 面積分送，驗證
全域覆蓋、唯一貢獻與 reference hashes。再做 HPC-07B 的 owned-cell traction
與 owner 膜更新，逐項核對合力、力矩及離散功；不能只比較平均壓力。
接著接上 HPC-07C 的全域 RMS／Aitken、共同失敗與 transaction 邊界。
最後 HPC-07D 用同一實際 FSI 問題比較單／多 rank，包含跨分區 triangle、空
rank、反向 traction、rollback／retry、checkpoint／restart 與收斂歷史。

原單分區 reference 保留，作為數值 oracle；移除拒絕條件只能發生在相應行為
與測試接線完成後。跨節點仍須實際 allocation，不能由單機 MPI 代替。
本報告只建立具體接線順序與驗收邊界，尚未宣稱任何分散式 FSI 求解已通過。

## 空 publication 分區已接線

`DistributedSurfaceLayout` 現在允許多分區中的 owned nodes／area vectors 同時為空；
全域 reference node count 仍須正、完整 reference positions／triangles 仍須合法，
單分區仍必須擁有全部 global nodes。非空 layout 的 hash stream 不變；空分區
也有由 partition rank 區分的 publication identity，不能交換 stamp。

契約測試驗證空 kinematics／traction、不同空分區 stamp 拒絕、非空 force 對空
layout 拒絕及 single-rank 空 layout 拒絕。MPI ownership 測試新增 rank 0／2
均空、rank 1 擁有三個節點的配置；覆蓋檢查通過，漏一節點共同拒絕，健康重試
通過。原 ownership 故障測試保持。三份 rank reports exit 0、無 timeout。

證據：`outputs/hpc07/empty-surface-v1/audit.json`。兩個測試目標均以 warnings
開啟重建；dependency-free contract 測試 exit 0。編譯與測試限制在 CPU 14／15，
MPI 三 ranks 共用這兩個 CPU，這是正確性測試，不是效能資料。

```bash
make -C solvers/coupling surface-contracts-test
make -C solvers/cpu parallel-ownership-test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
```

原表中的空 layout 限制已在本批解除；Aitken local weights 的正總和限制與各 FSI
單分區限制仍在。Reference 面積貢獻路由、ghost exchange 與真正分散式 fluid
traction 尚未完成，HPC-07A 仍未勾選。

## 空 rank 的 Aitken 局部貢獻

`DynamicWeightedAitkenRelaxation` 現在接受空 local weights 與空 iterate／residual。
空 inner product、numerator、denominator 均為 0，local residual scale 提供既有
reference floor；global weight total 仍須有限且嚴格為正，非空 weight 每項仍須
有限且正。空 rank 不得使用單 rank `Propose` 繞過 global reduction。既有非空
weight／state 的 hash stream 不變。

契約測試比較空 participant 與完整 owner 的 initial／dynamic proposal、pending／
accepted control identity、係數、reset；錯誤 residual size、foreign proposal、
global total=0、非空 zero weight 仍拒絕。原加權／不等分區與 stale proposal
測試保持通過。

三 rank ownership 測試另接真正 MPI SUM／MAX：rank 1 有兩個殘差分量，rank 0／2
為空，兩次 proposal 的 relaxation 與單 rank reference 完全相等，且全部 ranks
的 pending／accepted control identity 一致。三份 rank reports exit 0、無 timeout。
此處 reduction 位於測試 driver，不冒充正式 FSI coordinator 已接線；原強耦合
與膜的單分區限制仍在。

證據：`outputs/hpc07/empty-aitken-v1/audit.json`。重現沿用前節兩個 make test
目標。本次 compile／test 限於 CPU 14／15，背景活動也記錄於正在執行的 rank
scaling 目錄，不能因 CPU binding 分離就宣稱系統資源完全無干擾。

## Sparse owner 加總元件

`OwnedScalarContributions.hpp` 新增 `SumOwnedScalarContributions`。每 rank 宣告
owned UInt64 IDs 與 `(node_id, scalar)` contributions，按 ID 分送到 routing shard；
shard 檢查唯一 owner、累加後只把結果送到 owner，輸出沿輸入 owned ID 順序。
支援完整 UInt64、空 root、全空 communicator，以及 owned node 無貢獻時的 0。

不建立 global_nodes 長度的場，也不 gather 完整 owner catalog。沿用 bounded
Alltoallv wire exchange；預設每次每方向 64 MiB 與每 rank owner／contribution
合計一百萬 records。這是 wire／record cap，不是 RSS cap；ID 分布不均仍可
集中 shard。局部錯誤以 collective stages 傳播，未提供 MPI process-loss recovery。

Repeated contributions 是正常加總，不能由本元件判定是否重複物理積分。節點
完整覆蓋仍由 publication ownership validator 證明；三角形 authority／reference
面積組裝是下一步。加總用 long double，輸出 double 前先檢查範圍；不承諾跨
rank 配置逐位元相同。

三 rank 測試有空 root、rank 1 為唯一 owner、超過 Int64 的 ID、重複正負貢獻、
零貢獻 owned node 及全空輸入。六項負向條件（重複 owner、未知節點、非有限
貢獻、wire cap、record cap、sum 超出 double 範圍）均共同拒絕，逐項健康重試
通過。原 surface／Aitken／ownership 測試也通過。證據在
`outputs/hpc07/scalar-contributions-v1/audit.json`；最終三份 rank reports exit 0、
無 timeout。測試沿用 `make -C solvers/cpu parallel-ownership-test`。

## Reference 三角形面積已分散組裝

`DistributedSurfaceAreas.hpp` 的 `ComputeDistributedSurfaceAreas` 接收已驗證
reference layout 與本 rank 的 triangle indices。先驗證 node publication 全域
覆蓋，再驗證 reference triangle array 的每個 index 恰有一個 owner；triangle
owner 不必擁有任何相關節點。每個 triangle 計算 A/3 的三個 nodal contributions，
以 sparse owner 加總路由，僅返回本 rank 的 owned lumped areas 與 global area。

舊 layout 中的 area 數值不參與計算，輸入 layout／stamp 不修改。Caller 安裝新
weights 後須重新取得 partition identity；不得沿用依舊 weights 產生的 stamp。
目前輸入仍是合法、完整 replicated reference layout，尚未改成 ghost-only 幾何。
此處的 index coverage 不是 full manifold／self-intersection 認證，其仍由 reference
mesh 的上游建立流程負責。

每個 owned node 最終 area 必須正且有限；全域 node／triangle 面積以 long-double
SUM 比較，門檻為 128×double epsilon×global triangle area。只對兩個標量做
Allreduce，沒有全域 nodal area vector。Record cap 包含 owned nodes 加三倍本地
triangle 數，wire cap 與既有 sparse sum 相同；上游 replicated-reference ownership
validator 的 metadata 成本仍存在。

三 rank 測試涵蓋 rank 2 擁有 triangle、rank 0 空、rank 1 擁有全部 nodes；另有
分散 node ownership。非等面積兩 triangle 的面積分別 0.5／1.0，共享 nodes 的
權重均為 0.5，全域面積為 1.5。測試刻意給舊 layout 不同權重，確認使用實際
幾何重算。重複 triangle、漏 triangle 與 record cap 三項共同拒絕／健康重試均
通過。三份 rank reports exit 0、無 timeout，原 ownership／Aitken／scalar tests
亦通過；來源與 logs 在 `outputs/hpc07/reference-areas-v1/audit.json`。

重現沿用 `make -C solvers/cpu parallel-ownership-test`。這完成 reference area
貢獻路由元件；ghost 場交換、fluid traction integration 與正式 FSI runtime 接線
仍待完成，HPC-07A 尚未勾選。

## Owner-to-ghost tuple 交換元件

`OwnedPointValues.hpp` 新增 `FetchOwnedPointValues`，接受 owned UInt64 IDs、
interleaved tuples、requested IDs 與共同 component 數。Owner publication 與
查詢按 ID 路由至 shard，再回傳到 requester；不建立完整全域 owner map 或
nodal field。允許重複查詢並保留輸入順序，支援空 owner／requester 與全空群組。
每次呼叫重新發布輸入，沒有隱含快取。

Duplicate owners、無 owner 的 query、非有限 tuple、component 分歧與 cap
超限共同拒絕。沿用每方向 wire cap 與 local owner＋query record cap；它們不
代表 RSS 上限，ID skew 仍可能集中 shard。所有 owned tuples 都會發布，即使
本次沒人查詢；這項成本需納入之後實際 ghost 通訊量測。

三 rank 測試以空 root 查詢 rank 1 的場，rank 2 為空 requester，包含 UInt64
最大 ID、ID 0、重複且反向排序查詢及 signed zero，tuple bits 保留。全空群組
與六項共同錯誤／逐項健康重試均通過，原面積／scalar／Aitken／ownership 測試
保持。三份 rank reports exit 0、無 timeout，證據在
`outputs/hpc07/owned-point-values-v1/audit.json`。

本元件只處理 ID 與 tuple；expected material catalog 的完整覆蓋、FSI time／
step／coupling iteration／producer-state stamp 綁定仍由下一層驗證。不能用
本測試宣稱已拒絕 stale FSI kinematics；正式 ghost publication wrapper 尚待接入。
重現仍用 `make -C solvers/cpu parallel-ownership-test`。

## FSI kinematics stamp 與 ghost 交換

`SurfaceGhostKinematics.hpp` 的 `FetchSurfaceGhostKinematics` 在 tuple 交換前
驗證 publication 與 caller 提供的 expected interface／stamp 完全相符，包括
producer-state 與 partition identity。全群組另核對 interface、reference／layout、
time、step、coupling iteration；producer-state 可以是 partition-local identity，
不錯要求每個 owner 的局部 state hash 相同。空 owner／requester 也要通過這些
檢查。之後驗證 node ownership 完整覆蓋，再交換 displacement／velocity 六個
component，按 requested IDs 回傳。

Expected stamp 必須由 caller 的 transaction authority 獨立取得，不能從待驗證
publication 自行複製就宣稱拒絕 stale state。此 wrapper 不從場值推導或證明完整
producer solver state；它核對的是既有出版契約的身分。回傳為獨立 requested-node
資料型別，不能直接當成 aligned-owned `SurfaceKinematics` 發布；其使用生命週期
仍由外層 trial／epoch 管理。

三 rank 測試包含空 root 查詢、local owner 查詢與空 requester；不同 owners 的
producer-state hashes 合法且不同。七種情境均共同拒絕並健康重試：過期時間、
錯誤 producer state、錯誤 subsystem、錯誤 partition、單 rank 的 expected 與
publication 同時改 iteration（由全域 epoch agreement 抓到）、未知 requested
node，以及非有限 displacement。原 tuple／area／Aitken／ownership 測試通過。
三份 rank reports exit 0、無 timeout，證據為
`outputs/hpc07/surface-ghost-v1/audit.json`，重現沿用 parallel-ownership-test。

目前是已驗證的 stamped ghost exchange 元件；尚未移除正式 FSI runtime／patch
map 的單分區限制，也沒有將 traction 或膜求解分散。HPC-07A–D 保持原未完成範圍。


## 分區 reference patch map

`MaterialSurfacePatchMap` 現在允許空或部分 owned publication slice，各 rank
仍保存完整 reference positions、triangles 與 material vertex mapping。全域
reference identity 與映射不隨 ownership 改變，map identity 仍包含 partition
identity，因此不同分區的 map 不會混用。跨 rank 的唯一 ownership 與完整覆蓋
由 MPI caller 驗證；建立局部 map 本身不宣稱完成這項集體驗證。

`MaterialSurfacePatchKinematics::ComposeTarget` 明確保留完整單分區要求，避免
將部分位移套用至完整封閉曲面而靜默遺漏其他 owner。正式 runtime 的既有單分區
檢查也保持；下一步需接入經 stamp 驗證的遠端節點值，再處理完整 trial。

契約測試建立三個 partition maps（含空 root），核對所有 reference vertex／
triangle mappings 與單分區一致、reference identity 相同、partition map identities
相異。有效的部分 publication 先通過其 layout 驗證，再確認 serial composition
以指定單分區錯誤拒絕。原有 patch geometry／kinematics 契約測試亦通過。

重現：`make -C solvers/cpu material_surface_patch_kinematics_test` 後執行
`solvers/cpu/material_surface_patch_kinematics_test`。本次建置及執行綁定 CPU 14、15，
均 exit 0，編譯無新增 warning；來源、binary 與 logs hashes 記錄於
`outputs/hpc07/partition-patch-map-v1/audit.json`。這是 reference map 層的支援，
尚非正式分散式 FSI runtime，HPC-07A 保持未勾選。


## Triangle owner 的變形座標與速度

`DistributedSurfaceTriangleKinematics.hpp` 將 reference triangle ownership 接入
stamped ghost 交換。每個 rank 只查詢其 triangles 所需的不重複 nodes，並按
原 triangle 順序與 winding 回傳 reference index、node IDs、reference position
加 displacement 的當前座標，以及 nodal velocity。Triangle owner 可以完全
不持有節點；先共同驗證 triangle coverage，再沿用 publication ownership、
expected stamp 與全域 epoch 檢查。Reference geometry 目前仍完整複製。

此元件提供後續 P1 牽引力積分的幾何輸入，不驗證封閉曲面、自交、變形 Jacobian、
clamped seam 或 backward-Euler 一致性；這些仍是外層 material trial 的責任。
回傳資料的生命週期由 caller 的 transaction 管理，不能跨 epoch 隱含快取。

三 rank 測試以 rank 1 持有全部 nodes、rank 2 持有 triangle、空 root 不持有
任何資料，解析比較三個非均勻位移／速度的結果。重複 triangle、遺漏 triangle、
record cap 與空 root 的 stale publication 四項共同拒絕及逐項健康重試通過，
原 ghost／area／Aitken／ownership 測試也通過。三份 rank reports 均 exit 0、
無 timeout，CPU affinity 為 14、15；來源、binary 與 logs hashes 在
`outputs/hpc07/triangle-kinematics-v1/audit.json`。重現沿用
`make -C solvers/cpu parallel-ownership-test`。正式流體牽引力與 FSI trial 接線
尚未完成，HPC-07A/B 維持未勾選。


## 已積分角點力送回 node owner

`AssembleOwnedSurfaceForces` 接受唯一 triangle owner 已積分的三個 P1 corner
forces，依 reference triangle node IDs 加總至唯一 node owner。先驗證 node
ownership 與 triangle coverage，然後使用 `SumOwnedVectorContributions`；不將
component 編碼進 ID，因此 UInt64 最大 node ID 仍合法。Vector helper 暫用
三次 scalar routing，wire cap 適用每次交換，不能視為三次總 traffic 或 RSS
上限。所有輸入保持不變。

測試使用兩個共享邊的 triangles、非均勻三分量 corner forces。Rank 0／2 各
持有一個 triangle，rank 1 持有全部 nodes，逐點解析值完全相符。MPI 全域
合力 `(8,20,12)` N、原點力矩 `(6,-9,-2)` N m，以及測試速度場 `v=(x,y,1)`
下的離散功率 30 W 均精確相符。重複／遺漏 triangle、陣列長度錯誤、非有限
第三分量、record cap 與第三分量累加超出 double 範圍六項共同拒絕／重試通過。
既有 ownership、area、ghost、triangle kinematics、Aitken 測試也通過。

三份 rank reports exit 0、無 timeout，綁定 CPU 14、15，證據在
`outputs/hpc07/surface-forces-v1/audit.json`；重現沿用
`make -C solvers/cpu parallel-ownership-test`。此元件不從流體場積分 corner force，
也不以 lumped area 除法替代 consistent traction mass projection。Caller 仍需
將角點力綁定到正確 trial／材料幾何，並將跨 cut-cell 的積分先送到 triangle
owner。這些正式接線、投影與結構更新仍待完成，HPC-07B 只記部分進度。


## 依 fluid cell authority 直接分送力

`AssembleOwnedSurfaceCellForces` 接受完整 dense background cell catalog 的唯一
ownership，每個 owned cell 提供已積分的 `(material node ID, force)` 列表；
沒有 retained interface quadrature 的 cell 也須以空列表參與 coverage。這讓
同一 material node 或 triangle 跨多個 fluid cells 的所有貢獻直接加總到 node
owner，不必先集中到 triangle owner，也不會按 triangle ID 去重。

元件核對全域 cell count、唯一 cell coverage、surface node publication ownership
與 vector contribution 有限性。它不證明 cell 內 quadrature 完整性、catalog
provenance 或 trial stamp；這些必須由後續 fluid integration caller 綁定。
Local cell count 與 owned nodes 加 nodal contributions 各受 record cap 限制；
vector sum 仍走三次 scalar routing。既有 triangle-corner API 保留給已集中於
triangle owner 的資料，正式 owned-fluid-cell 接線可直接使用新 API。

三 rank fixture 使用兩個有力貢獻的 cells 與一個空 contribution cell；兩個
cells 的共享 node forces 和既有解析 nodal values 相同，原合力／力矩／功率
檢查保持通過。重複、遺漏及越界 cell、未知 material node 與非有限第三分量
五種錯誤共同拒絕並健康重試。其餘 parallel ownership 回歸通過，三份 reports
exit 0、無 timeout，CPU affinity 14、15；證據在
`outputs/hpc07/surface-cell-forces-v1/audit.json`。重現沿用
`make -C solvers/cpu parallel-ownership-test`。HPC-07B 尚未完成正式流體積分、
consistent projection 與結構更新驗收。


## 共用 P1 force／consistent mass 單點核心

`SurfaceP1TractionContribution.hpp` 將單一 quadrature point 的
`traction_component * Ni * weight` 與 `Ni * Nj * weight` 整理為三個 corner
forces 和 3×3 consistent mass block。現有 `BuildFluidSurfaceTraction` 已使用
此核心，原 Cauchy stress、canonical corner provenance、乘法順序、force
compensated accumulation、mass 累加及 projection identity stream 保持。

此核心驗證輸入與乘積有限性，不自行限制 signed quadrature weight 或重新
正規化 barycentric values；幾何與 shape provenance 仍由 catalog authority
驗證。後續分散式 cell integration 可共用相同公式，同時保留 force 與 mass，
不改成 lumped projection。

三點 degree-two triangle rule 對線性 nodal traction 的解析驗證通過：面積
0.5 的 consistent mass 對角 1/12、非對角 1/24，nodal force 等於解析 mass
乘 nodal traction。非有限 weight／shape／traction、force overflow 與 mass
overflow 五項拒絕通過。原 fluid traction 壓力、黏性、非均勻場、抵銷、平移
及 publication authority 回歸全部 exit 0，無 timeout，編譯未啟用 NDEBUG。
證據為 `outputs/hpc07/p1-traction-kernel-v1/audit.json`，重現命令為
`make -C solvers/cpu fluid-surface-traction-test`（指定本機 PETSC_DIR）。
這次驗證針對單分區正式 traction consumer；分散式 cell integration 與
consistent mass projection 接線仍未完成。


## 單 owner 的 bounded consistent projection

`ProjectDistributedSurfaceTraction` 接受各 rank 的 `(row material ID, column
material ID, mass entry)` 與已加總的 owned nodal forces。Mass entries 路由至
指定 projection owner，force tuples 亦按 node IDs 取到該 owner；只有它配置
完整 dense mass matrix，使用既有 consistent Cholesky solver，然後將 traction
送回各 node owner。Projection owner 與 node owner 可以不同，不要求 root
持有 material nodes。這符合先集中小型 projection 的契約，尚非分散式線性求解。

預設最多 4096 nodes，dense double matrix 本身最多 128 MiB；wire cap 限制
每次交換，不能解讀成包含 packets、layout、matrix 與 solver vectors 的 RSS
上限。Replicated reference geometry 仍在。核對 node coverage、projection owner
agreement、dimensions、有限 mass／force、未知 IDs、mass symmetry 與 Cholesky
正定性；反覆同一 matrix entry 是合法的跨 cell 積分加總。Caller 仍須證明
force／mass 來自同一 trial 與唯一 quadrature authority，代數元件不自行產生
可信的 publication stamp。

三 rank 製造解使用共享邊的兩個 P1 triangles，rank 0／2 分別提供局部 consistent
mass，rank 1 持有全部 nodal RHS。解析 traction 為 `(1+x,2+y,3+x+y)`，分別
選 rank 0、2 求解，全部 component 誤差小於 1e-13，空 node owners 回傳空場。
零 mass、未知 mass node、非對稱 matrix、node cap 與 owner 分歧五項共同拒絕
並健康重試；完整 ownership 回歸也通過。三份 reports exit 0、無 timeout，
整套測試約 15.6 秒，RSS 分別 38,879,232／38,739,968／38,662,144 bytes。
這些是整套正確性測試數據，不能當成單獨 projection 效能或大型記憶體證據。

來源、binary 與 logs hashes 在 `outputs/hpc07/surface-projection-v1/audit.json`，
重現沿用 `make -C solvers/cpu parallel-ownership-test`。正式 fluid quadrature、
projection stamp、膜更新及 FSI runtime 接線仍待完成，HPC-07B 尚未勾選。


## 同源 P1 force／mass 的分散式組裝

`AssembleDistributedSurfaceTraction` 現在由各 owned cell 的 P1 point records 同時
建立 corner forces 與 consistent mass entries，再核對完整 fluid cell ownership，
將 force 加總到 node owners，經 bounded projection 回傳 owned traction。每個
point 的 node IDs、barycentric values、traction 與 weight 同時供兩種積分使用。
零 interface points 的 cell 仍參與 coverage，mass records 在配置前檢查上限。

三 rank 測試以兩個共享邊 triangles 的 degree-two 積分及線性 traction 製造解
串接整條代數路徑。Owned RHS 對解析值誤差小於 1e-14，投影 traction 小於
1e-13；無效 weight、重複 cell 與 mass record cap 三項共同拒絕／健康重試通過，
其餘 ownership／projection 回歸保持。三份 reports exit 0、無 timeout，來源及
logs hashes 位於 `outputs/hpc07/traction-assembly-v1/audit.json`，重現沿用
`make -C solvers/cpu parallel-ownership-test`。

此輸入仍為 caller 提供的 stress／provenance records；下一步必須從實際 IGA
velocity／pressure、surface quadrature catalog 與 trial authority 建立它們。
回傳刻意是無 stamp 的 owned values，不能冒充已綁定材料狀態的 SurfaceTraction
publication。正式 FSI runtime 與 HPC-07B 驗收仍待完成。


## 實際 owned-cell IGA stress 與 catalog 擷取

`BuildOwnedFluidSurfaceTractionPoints` 由 owned fluid cells、局部 retained-cell
IGA coefficients、material／domain／surface catalog 與 patch map 建立 P1 records。
它在實際 quadrature parametric positions 評估 IGA pressure／velocity gradient，
呼叫既有 `FluidOnStructureCauchyTraction`，保留 canonical barycentric corner
順序，按 source vertex provenance 對應 material node IDs。只積分本地選定
patch 的 points；沒有 retained patch points 的 owned cell 回傳空列表。

驗證 material／topology／domain／catalog identity、viscosity、cell 範圍與重複、
retained state 完整性、係數大小及有限性、point label／source triangle／node
mapping；多餘 state 也拒絕。這是 local extractor，MPI caller 必須以 collective
local stage 協調錯誤，並驗證跨 rank 共享 IGA coefficients 和 trial authority。
目前仍保留 replicated material／domain／catalog；未宣稱幾何已完整分散。

1／3 ranks 的實際 tetrahedron patch 案例，按 cell ID 分配積分工作與 material
node 分配發布權，串接 P1 assembly／bounded projection。壓力案例的 traction
為 `(0,0,-20)` Pa；affine viscous 案例為 `(6,0,8)` Pa。各 owned force／traction
與既有 serial oracle 在原尺度下 1e-12 容差內相符。兩案例均演練缺少 retained
cell state、非有限 pressure coefficient、係數數量錯誤，所有 ranks 共同拒絕，
健康重試後仍通過數值比較。原 fluid traction 全套回歸也通過。

最後四份 rank reports 均 exit 0、無 timeout，來源、binary 與 logs hashes 在
`outputs/hpc07/owned-iga-traction-v1/audit.json`。首次 launcher 因缺少 report 父
目錄而未啟動測試程式，補建後的初版與最終版測試均成功，此錯誤列於 audit。
重現：建置 `fluid_surface_traction_test` 後以 `mpiexec -np 1` 與 `-np 3` 執行。
測試為建立 oracle 會在每 rank 計算完整 serial reference；新 extractor 僅對
owned cells 積分，不能把測試總耗時當成分散式求解效能。

這批接入了實際 IGA stress 與 catalog，但正式 runtime 的 epoch／producer-state
stamp、共享 coefficient 一致性、守恆診斷與膜更新仍待完成。HPC-07B 未勾選。
