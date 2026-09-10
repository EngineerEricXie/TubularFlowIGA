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
