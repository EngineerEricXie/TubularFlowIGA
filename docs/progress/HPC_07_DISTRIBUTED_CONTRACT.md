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
