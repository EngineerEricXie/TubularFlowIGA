# HPC-03D：移動幾何與 distributed extension 進度

HPC-03D 尚未完成。正式 moving flow 仍使用 PETSC_COMM_SELF 與完整 global
state，固定幾何 ImmersedTransientDistributedRuntime 尚不能跨 active layout
凍結 history。此處追蹤以固定 Cartesian background、穩定 node IDs 接入真正
分散式 history extension；不把每個 rank 複製完整求解視為完成。

## 共用幾何拓樸建立

ImmersedVelocityExtension 現在把 InitializeTopology 與 dense matrix／factor、
velocity／pressure extension 分開。拓樸階段只讀 old／new geometry、active
layouts、extension layers 與 limits，保留 positive-cell／band／forward-reverse
coverage、anchor／unknown IDs、face traces、anchoring checks；不讀 numerical
state，也不配置 dense matrix。數值 state identity 與既有 serial dense solve
仍由原 Build 路徑負責。這是後續 MPI sparse operator 共用原方程的準備，
尚未實作分散式 extension。

拆分前後均重建並執行 immersed_velocity_extension_test。測試新增輸出四個
identity 作為比較依據；兩次完整 stdout byte-identical：

- Extension: `f0ce30930dca32264a572600e6c206a8c18975cd353f776fa6bf0a67714388ec`
- Operator: `3dca1863bd4133b4e0a96753ae982775c1b3c065cab1de575804f780cd1311b4`
- Reduced matrix: `0a057e7377a7526bf47df7bff76365f4b8ba6461495572f334ff17550c565399`
- History: `8f95890a44efe635e7177eec190e3dece3a45743121ed150fda5fd31224a291c`

報告 fixture 有 88 band cells、517 nodes、345 unknowns、184 faces；polynomial
relative max 為 1.3777867735598193e-13，獨立 LDLT relative L2 為
5.8578413540129128e-11，normalized residual max 為 3.2212246812575427e-20。
既有 stationary／contraction、cap、fault、anchor bitwise 與 retry assertions
也全部通過。四個輸出 identity 僅代表該報告 fixture，其他案例以原測試判斷。

第一次 baseline build 暴露測試缺少明確 PrescribedSurfaceMotion.hpp include；
已修正 include 與 Makefile dependency，再完成 baseline 與 refactor build/run。
原失敗 log 保留。證據在 outputs/hpc03/moving-extension-topology-v1/audit.json，
同目錄保留原 header、baseline hashes、完整 build logs 與相同的結果 logs。
兩次成功 build 無 compiler warnings；可用下列指令重跑現況：

```sh
make -C solvers/cpu immersed-velocity-extension-test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
```

後續須完成 owned sparse extension assembly／solve、必要 anchor 與 unknown
halo exchange、target history ownership 更新，並接入 moving runtime 的 prepare／
commit／rollback。幾何與拓樸仍複製的記憶體成本需另外量測與追蹤。


## MPI sparse extension

DistributedImmersedVelocityExtension 已共用上述拓樸與 third-normal-derivative
face traces，以穩定 node ID 建立 unknown rows。每 rank 只組裝其 owned sparse
rows，anchor tuples 透過 sharded lookup 取得；所有 source anchors 的唯一
ownership／coverage 也會驗證。原始與 diagonal-scaled matrix 都是 PETSc
分散式矩陣，四個 RHS（velocity xyz、pressure warm start）使用 MUMPS LU，
可透過 private `immersed_extension_` options snapshot 配置。

求解後回到原始矩陣檢查
`||A x − b||₂ / (||A||∞ ||x||₂ + ||b||₂)`，必須小於 extension residual tolerance。
只傳回 caller 要求的 target tuples，committed anchors 直接沿用原始 double bits。
沒有 unknown 的 stationary 路徑不建立矩陣或 solver，仍驗證 source coverage。
這不是 pressure 方程；第四欄僅為後續求解的 scalar warm start。

6×6×6 Cartesian fixture 的 tetrahedron 從 x=0.18 移到 x=0.40，active node
set 確實改變。串行 dense oracle 僅在測試 rank 0 建立，分散式 extension 本身
不建立完整 numerical state 或 dense factor。最終測試結果如下：

| ranks | relative L2 對串行 | 最大絕對差 | retry relative max | 原始系統 normalized residual max |
|---|---:|---:|---:|---:|
| 1 | 3.24313e-14 | 5.93303e-13 | 0 | 1.96804e-17 |
| 2 | 3.21620e-14 | 6.28830e-13 | 6.56850e-16 | 2.65551e-17 |
| 4 | 3.34372e-14 | 6.28830e-13 | 1.80165e-15 | 2.34095e-17 |

Source anchors 在求解與 retry 均 bitwise 相同，包含 signed zero。缺少 source、
nonfinite coefficient、非法 target ID、Close 後使用皆共同拒絕。另一個時間
推進但幾何不變、layers=0 的案例驗證零 unknown；全部 source 位於 rank 0，
其餘 1／3 ranks 為空 source owners，分散 query 結果與原始 fields bitwise 相同。

最初雙 rank 測試對重複 MUMPS solve 要求 bitwise 相同而失敗；已保留失敗 logs，
改為記錄並要求 retry relative max < 1e-12，同時保留嚴格 anchor bitwise gate。
最終 1／2／4 ranks 共七份 reports exit 0、無 timeout，編譯無 compiler warnings。
證據在 `outputs/hpc03/distributed-extension-v1/audit.json`，可重跑：

```sh
make -C solvers/cpu distributed-immersed-extension-test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
```

限制：geometry／face traces 仍 replicated，owned sparse maps 也暫時與 PETSc
storage 並存；尚未量測大型案例的記憶體收益或 scaling。此 class 產生數值 tuples，
committed epoch authority 由 caller 提供；PETSc Vec 抽取、target history provenance／
clock、移動 runtime 的 prepare／rollback 接線仍未完成。HPC-03D 保持未勾選。
