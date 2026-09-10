# HPC-03D：移動幾何與 distributed extension 進度

HPC-03D 尚未完成。分散式 sparse extension、mapped history、moving-wall
operator 與擁有 committed／trial epochs 的 MPI Newton runtime 已接線。
連續兩步改變 active nodes 的 1／2／4-rank runtime 串行場對照、逐位元
提交、失敗重試與跨 epoch options snapshot 檢查已通過。正式 graph／FSI
moving flow 尚未改用這個分散式 runtime。

Production numerical vectors／matrices 維持分散；geometry／face traces 目前
仍 replicated。此報告的測試 oracle 可建立小型完整串行場，不屬於 production
計算路徑。下列 audits 各自綁定當時的 source／binary，不代表同一個 HEAD，
後續修改前保存的版本用於重核 hashes。所有結果都是單工作站 MPI 驗證，
未宣稱大型案例記憶體收益、scaling 或跨節點能力。

## 共用幾何拓樸建立

ImmersedVelocityExtension 現在把 InitializeTopology 與 dense matrix／factor、
velocity／pressure extension 分開。拓樸階段只讀 old／new geometry、active
layouts、extension layers 與 limits，保留 positive-cell／band／forward-reverse
coverage、anchor／unknown IDs、face traces、anchoring checks；不讀 numerical
state，也不配置 dense matrix。數值 state identity 與既有 serial dense solve
仍由原 Build 路徑負責；以下 MPI sparse extension 共用此拓樸。

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

幾何與拓樸仍複製的記憶體成本需另外量測與追蹤。


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
committed epoch authority 由 caller 提供。下述 runtime 使用其自身 committed
Vec 與 accepted clock 呼叫 extension，沒有把傳入 clock 當作獨立權威。


## 從 committed PETSc Vec 讀取 extension source

ExtendCommitted 直接讀取來源 Vec 的 owned node rows，不建立 gather 或全域
numerical state。它核對 constructor 綁定的 source active-layout hash、Vec size／
communicator、完整四欄 node ownership，以及 source／target geometry time 與
共同 step index／dt。Clock 需滿足 target_index=source_index+1 及精確的時間
推進關係。來源 Vec 是否確為 runtime committed state，仍由呼叫它的 runtime
負責證明，傳入 clock 本身不提供權威身分。

抽取只包含 velocity xyz 與 pressure warm start；controller／gauge rows 不納入。
Owned coefficients 複製至 local staging，再使用既有 sparse extension 與 target
query 路徑。來源 Vec 從未寫入，失敗時也不需回復其內容。

1／2／4 ranks 以 contiguous whole-node Vec ownership 驗證，與前一個 round-robin
source tuple 測試採不同分配。對 serial dense oracle 的 relative L2 分別為
3.24313e-14、3.21605e-14、3.28674e-14；所有 source Vec coefficients 與 target
anchors 逐位元保持，包含 signed zero。Null source、錯誤 source layout／time
共同拒絕；多 rank 另驗證各自有效但不一致的 clock、PETSC_COMM_SELF 來源，
以及切斷四欄 node 的 Vec ownership。拒絕後健康重試成功。

另以 stationary／zero-unknown 案例，把全部 node rows 放 rank 0，controller
與 gauge rows 放最後一個 rank。多 rank 的最後一個 owner 沒有任何 node rows，
仍能正常參與 extension；兩個 scalar rows 原值保持，沒有混入 nodal output。
此案例驗證排除 scalar rows；controller 搬移由下述 target-seed 案例驗證。

最終七份 rank reports exit 0、無 timeout，build 無 compiler warnings。證據在
`outputs/hpc03/committed-extension-v1/audit.json`，重現沿用
`make -C solvers/cpu distributed-immersed-extension-test PETSC_DIR=...`。



## Target history 與暫態體積組裝

ImmersedDistributedVelocityHistory::FreezeMapped 已透過 extension 的精確 target
layout／communicator 與 committed-source clock 檢查，取得 target cells 所需的
halo。History 僅保存三個 velocity fields，並保留 source geometry identity、
source／target time/index，以及每個 required node 的 Committed／Extended
provenance。候選 fields 與 metadata 全部通過共同驗證後才交換發布；既有固定
layout 的 Freeze 仍使用原 VecScatter 路徑，並標記全部為 Committed。

ImmersedDistributedTransientVolume::Freeze 可選擇 mapped source，先共同核對
各 rank 使用同一種模式，再凍結 history 與 body force。Force callback 失敗會
撤銷新 history；對 active inputs 再次 Freeze 則拒絕並保留原 inputs。此處只接入
單一 immutable target geometry 的體積積分。完整固定幾何 operator 的零壁速
契約保持；moving-wall 方程由獨立 moving operator 接入。

1／2／4 ranks 的移動 active-set fixture 已以 mapped history 實際組裝 Jacobian
及 residual，對串行 dense extension 所提供的 local history，最大絕對差均為
6.77626e-21。Provenance 與 source／target clock 正確；freeze 後修改 source Vec
不會改變已凍結 halo 或組裝結果。由擁有實際 quadrature points 的 rank 注入
body-force 失敗，所有 ranks 共同清理，解除失敗後 retry 通過。

新路徑七份與既有 fixed history／volume 各 11 份 rank reports 全部通過：
1／2／4 ranks 與各自四 rank split communicator，共 29 份 exit 0、無 timeout。
Fixed volume 同時覆蓋 expanded／compact quadrature。證據位於
outputs/hpc03/mapped-history-v1/audit.json；執行期間保留建置來源 manifest，
後續 target-seed 修改前另保存三份原 source／binary，所有 hashes 核對一致。
最初 test callback 編譯有一個 misleading-indentation warning，已在正式測試前
修正，final build 無 compiler warnings。



## Target state seed 與 controller 搬移

BuildOwnedTargetSeed 依 target Vec 的 ownership 回傳 local row-order candidate，
不寫入 source 或 target Vec。它要求精確 target layout／communicator、合法 Vec
size 與完整 node ownership，沿用 ExtendCommitted 的來源及 clock 檢查。四個
nodal fields 使用 sparse extension；controller 依穩定 UInt64 ID 從 source owner
取得，target gauge 設為 canonical positive zero，與 serial moving seed 政策一致。
Target 所需但 source 不提供的 controller 由 owner lookup 拒絕。

1／2／4 ranks 的移動 active-set seed 對串行四欄係數最大絕對差分別為
5.93303e-13、6.28830e-13、5.86198e-13。另一個 stationary fixture 將 controller
從最後一個 source rank 搬到 root target owner，其餘 ranks 的 target ownership
為空。Nodal values 逐位元保留、controller 值正確、gauge 為正零。非有限 source
controller 共同拒絕，修正後 retry 完全一致；target Vec 中預放的 sentinel 保持。

七份 rank reports 全部 exit 0、無 timeout，build 無 compiler warnings；同時重跑
mapped history／volume 與既有 extension 案例。證據為
outputs/hpc03/target-seed-v1/audit.json，沿用 distributed-immersed-extension-test
指令。此方法只建立候選 vector，clock／state 的提交權責留給 runtime。


## Moving material-wall operator

ImmersedMovingTransientDistributedOperator 已建立單一 immutable moving endpoint
的分散式組裝入口。它共用 fixed operator 的 sparse rows、volume、ghost、port、
gauge 及 Nitsche 方程；壁速由該 geometry epoch 的 canonical triangle／barycentric
provenance 取得。Constructor 共同核對 material content，Freeze 核對 material
step interval／evaluated time，並接入 mapped history。Fixed operator 仍使用零壁速。

組裝仍由 cell owner 執行，數值 state／matrix 維持分散；Eulerian field flux
診斷沿用共用 machinery，尚不代表已完成 moving geometric conservation 診斷。
此 operator 不擁有 accepted clock 或前一個 epoch，兩者由 moving runtime
管理。既有 all-flow port compatibility 規則保持，沒有擴張串行 API 的模式範圍。

目前 4×3×3 背景、平移 cube、material wall speed 0.16 的單 rank flow 模式已
與串行 BeginMovingTrial 完整組裝比較通過：residual relative L2 1.82272e-16、
Jacobian-action relative L2 4.38496e-16，physical diagnostic 最大差 3.38813e-21。
測試中的 serial oracle 用於比較，並非 production operator 的計算路徑。

第一次沿用 4×1×1 背景時，移動後沒有 full-cell ghost anchor，串行參考在
coverage gate 拒絕；擴充為 4×3×3 後保留相同 gate。第二次發現原測試預設
每個 cell 都含 port，新案例有 interior cells，已修正 oracle 只累加有該 label
的 quadrature。兩次失敗 logs 保留於 outputs/hpc03/moving-operator-v1。

五種模式 flow／pressure／mean-normal-traction／closed／inertial 的 1／2／4 ranks，
以及 fixed pressure 1 rank／inertial 4 ranks 回歸已全部通過：40 份 rank reports
exit 0、無 timeout。Across-suite residual relative L2 最大 3.85244e-16、
Jacobian-action relative L2 最大 2.54915e-15、physical diagnostic 最大差
1.08802e-13。196 份建置來源／binary hashes 與各 rank logs 核對通過，證據為
outputs/hpc03/moving-operator-v1/audit.json；後續修改前保留四份 tested 檔案。
上述 operator 比較不替代 runtime 的交易與守恆驗收。


## 固定背景 cell 分區

共用 sparse operator 新增 FixedBackground 模式，以完整背景 cell ID 的連續
區間分配 owner；inactive cells 也計入分區，邊界只由背景 cell 總數與
communicator size 決定。Volume／port／gauge／measurement 隨所在 cell，
ghost face 隨 canonical minus cell，純量 target 留在最後 rank。既有 CellCount
與 WeightedContiguous 模式保持。這是 assembly work 的固定分區；active
node rows 仍依當前 layout 分配，跨 epoch 數值由 extension／seed 路徑搬移。

新增測試將背景擴為 8×3×3，cube 平移 0.37（壁速 2.96），要求 source／target
active node IDs 確實不同，並核對兩個 epoch 的 owned cells 均位於同一背景
區間。完整組裝仍與串行 moving oracle 比較；非對齊案例單 rank 已通過，
residual relative L2 2.03897e-16、action relative L2 4.08256e-16，
physical diagnostic 最大差 7.11508e-20。Flow 的 1／2／4 ranks、其他四種
模式的 4 ranks，以及 4-rank split communicator 均通過，共 27 份 reports，
exit 0、無 timeout。最大 residual／action relative L2 為 3.5137e-16／
2.76832e-15，physical diagnostic 最大差 1.06803e-13。證據為
outputs/hpc03/fixed-background-v1/background-audit.json。

首次使用 0.4 平移時，串行 reference preflight 遇到 certified-empty cut-cell
quadrature 而拒絕。失敗 source／binary／logs 保留於 fixed-background-v1，
先以 0.37 的非對齊位置驗證固定 ownership，再追查此 preflight 失敗。

追查發現串行 ImmersedTransientFlowRuntime::PreflightCatalogs 已驗證空規則，
卻在判斷 point count 為零之前呼叫僅接受正體積規則的 iterator，令原本的
empty-cut-cell 分支無法執行。已將空 count 判斷移到 iterator 前，仍要求
ValidateUsableRule／ValidateUsableCompactRule 通過，且 Inside cell 為空仍拒絕。
原始 0.4 平移 fixture 以獨立 binary 重跑 1／4 ranks 全部通過：5 份 reports，
最大 residual／action relative L2 為 3.55513e-16／2.23435e-15，physical
diagnostic 最大差 1.95399e-14。證據為 fixed-background-v1/aligned-audit.json。
正式測試入口新增 aligned 選項以保留此案例，重新建置後單 rank 入口亦通過
（fixed-background-v1/cli-audit.json）；residual／action relative L2 為
1.71031e-16／3.95294e-16。


## Moving runtime epoch 交易接線

ImmersedMovingTransientDistributedRuntime 擁有獨立 committed／trial geometry
與 PETSc state，使用 FixedBackground assembly 分區。BeginTrial 核對 accepted
clock，從 runtime 自身 committed Vec 建立 target seed，凍結 mapped history／
force，再共同關閉 extension。求解期不依賴舊 epoch 或 extension 的數值資料；
失敗 Newton 可從 mapped seed 重試，AbortTrial 丟棄整個候選 epoch。

PrepareCommit 保留原 committed state；FinalizeCommit 僅交換已準備好的
state／epoch pointers 與 scalar clock，無 allocation、MPI 或 PETSc destroy。
舊 epoch 暫存為 retired，於下一次 collective BeginTrial 或 Close 釋放。

兩步求解、prepare failure、abort／retry、prepared rollback 的 1／2／4-rank
測試全部通過，共 7 份 reports，exit 0、無 timeout，195 份 source／binary
hashes 核對。證據為 outputs/hpc03/moving-runtime-v1/clock-audit.json。
這組平移案例未改變 active node IDs。第一次測試的第二步誤把 Evaluate 的
step_end 參數傳成 dt，修正為 start+dt 後重跑；原失敗 source／binary／logs
仍保留。

## Runtime 場對照與 exact publication

測試 root 建立 serial moving runtime，使用相同 geometry／flow／extension
設定，逐步比較速度三分量、壓力與 gauge 的 owned rows。每欄要求
`L2(error) / max(1, L2(reference)) < 1e-8`。小型全場 broadcast 僅用於測試。

每次提交後，prepared trial 的 owned field bytes 必須完整成為 committed Vec；
geometry／publication／layout identities 與 accepted clock 必須一起對應該
trial。另注入 rank-local body-force freeze failure、prepare failure，檢查
候選不發布、舊 committed Vec／clock 不變，以及 abort／retry。

未改變 active nodes 的 normal 模式已通過 2 ranks：兩步合計 Newton／KSP
各 4 次，最大分欄 L2 指標 7.93062e-15，exact publication 通過。
兩份 reports exit 0、無 timeout，195 份 source／binary hashes 與 logs 核對；
證據為 outputs/hpc03/moving-runtime-v1/publication-normal-audit.json。

changing 模式使用 8×3×3 背景，兩步各平移 0.37，每步都要求 active node IDs
改變。第一次對照在 serial oracle 的 KSP 失敗：serial 原設定為 FGMRES 且
未明確指定 LU backend，distributed 使用 GMRES／MUMPS。匹配後 KSP view
確認 gmres／mumps、relative tolerance 1e-10、absolute 1e-50，第一步場對照
通過。第二步被一層 extension forward-band gate 拒絕：0.37 位移大於
0.3 cell 寬，新增正體積 cell 超出該範圍。

兩邊都明確使用兩層 extension，另保留一層不足時必須拒絕且 committed
state／clock 不變的 negative case，再以兩層重試。早期單 rank 通過，
但雙 rank 第二步 1153-row 系統回報 KSP reason -11（PC_FAILED）。加入
PC failure reason 訊息後確認為 3（factor memory exhausted），原失敗 logs
保留於 moving-runtime-v1 的 band-*／pc-reason-2 目錄。

修復矩陣 options 綁定後，明確設 MUMPS ICNTL(14)=100 的 bound-runtime-test
通過 1／2／4 ranks，最大分欄 L2 指標分別為 2.71512e-13、4.05808e-13、
2.71825e-13；兩步合計 Newton／KSP 各 6 次。KSP view 核實 workspace 設定，
全部七份 reports 與來源 hashes 核對於 bound-audit.json。沒有加入 pivot
shift、修改矩陣或放寬 tolerances；workspace 裕量是這個案例的顯式設定。

## 矩陣選項與跨 epoch 快照

PetscSolverOptions::Attach(Mat) 將 private options database／prefix 綁定到矩陣，
讓 factor backend 讀取 runtime-prefixed 選項。共用 immersed Newton 在首次
組裝完成後、KSPSetOperators／factor setup 前綁定一次；serial transient
在 constructor 綁定已建立的 Jacobian，distributed extension 綁定 scaled matrix。
首次將 distributed 綁定放在 constructor 時被未組裝 guard 拒絕，修正綁定
時機後通過，沒有放寬 guard；workspace-2 失敗證據保留。

矩陣選項隔離測試通過 1／2／4 ranks 及內部 split communicators，共七份
rank reports。兩份同時存在的 private snapshots 分別指定 ICNTL(14)=37／81，
直接讀取 MUMPS factor 驗證，並驗證來源後續改動及 unprefixed 值 5 不會
覆蓋快照。原 options suite 同時通過，證據為
outputs/hpc04/matrix-options-v1/audit.json。其他 runtime 的 matrix-backend
選項尚未逐一稽核，這些結果不代表全專案矩陣選項驗收。

Moving runtime constructor 擁有 immutable options source；之後 Newton 與
extension epochs 都從它建立各自快照，使用旗標逐層傳回來源，epochs 先於
來源銷毀。共用 Newton 保留原 constructor，另提供顯式 borrowed source
入口，一致性檢查針對實際來源。第二步測試將各 rank 的 global moving／
extension KSP types 改成不同且無效的值，要求新 epoch 仍正常求解，測試
scope 結束恢復全域設定。

既有 fixed transient runtime 的雙 rank 回歸通過（field relative L2
1.11329e-14）；standalone extension 四 rank 回歸通過，並由 KSP view 驗證
其 private ICNTL(14)=73。共六份 reports 及 source hashes 核對於
moving-runtime-v1/compat-audit.json，覆蓋原 constructor 與獨立 extension API。

最終 runtime 驗證使用同一 binary 的 1／2／4 ranks，共七份 reports exit 0、
無 timeout；最大分欄 L2 指標分別為 2.71512e-13、4.05808e-13、2.03536e-13，
兩步合計 Newton／KSP 各 6 次。每次 active node IDs 均改變，freeze／prepare
failure、band rejection／retry、abort／rollback、exact publication 與 global
options 凍結檢查均通過。證據為 moving-runtime-v1/final-audit.json，核對
最終 source／binary 與各 rank logs hashes。雙 rank 的 snapshot-runtime-test
與最終重建 binary 的 SHA256 相同；其間只有兩個新 header 的空白格式調整，
另由 format-audit.json 保留核對。此驗證不替代 moving conservation 驗收。

## 重現目前的 runtime 測試

在工作站或已分配的 compute node 上執行；需要 real64 PETSc 與 MUMPS。
下列 runtime command 明確使串行 oracle 使用與分散式 core 相同的
GMRES／MUMPS，並明確增加 factor workspace 裕量；沒有變更 nonlinear 或
linear tolerances。

```sh
make -C solvers/cpu immersed_moving_distributed_runtime_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
for ranks in 1 2 4; do
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 timeout --kill-after=5s 900s \
    mpiexec --bind-to core --oversubscribe -np "$ranks" \
    solvers/cpu/immersed_moving_distributed_runtime_test changing \
    -moving_oracle_pc_factor_mat_solver_type mumps \
    -moving_oracle_ksp_type gmres -moving_oracle_ksp_converged_reason \
    -immersed_moving_mat_mumps_icntl_14 100 \
    -moving_oracle_mat_mumps_icntl_14 100 || exit "$?"
done
```

Production API 的 constructor 接收並擁有 initial geometry。BeginTrial 接收
新的 geometry ownership、target index 與 extension layers；失敗也會消耗該
候選 geometry，重試時需重新建立候選。所有 mutating operations 都需所有
communicator members 共同呼叫；FinalizeCommit 是 collective preparation
之後各 rank 的本地非拋例外發布。Close 與 runtime destruction 必須在
PETSc／MPI finalize 前完成。Caller 應依位移範圍選擇 extension layers，
超出範圍時會拒絕，不會默默縮減 target active set。

## 剩餘驗收

仍需完成 moving material／geometric conservation 與相應 failure gates，
接上正式 graph／FSI adapter 及跨 domain paired commit。Restart、大型案例的
記憶體／負載量測與是否需要動態重分區亦未完成，分別連同 HPC-05／06／07／09
追蹤。HPC-03D 仍未勾選。
