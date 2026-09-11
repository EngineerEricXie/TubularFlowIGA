# HPC-03D：移動幾何與 distributed extension 進度

HPC-03D 尚未完成。分散式 sparse extension、mapped history、moving-wall
operator、MPI Newton runtime、graph adapter 與 prescribed-motion case factory
已接線。兩步 active-set 改變的 runtime／graph adapter 已通過 1／2／4 ranks
串行場對照、提交隔離、失敗重試與 options snapshot 檢查。正式 CLI 的
零場啟動已通過 1／2／4 ranks、dt=.0625 的四步 graph／守恆驗收與端口
數值比較；dt=.125 的第二步不收斂紀錄保留。

薄切割 seed 支撐範圍、NNLS dual stopping、fitted quadrature ordering 與
fraction diagnostics 修正均有獨立證據。案例產生與保存結果比較工具已
交付，命令見 [MOVING_IMMERSED_CASE.md](../MOVING_IMMERSED_CASE.md)。
FSI owned material composition／traction publication 已有元件與 runtime
驗證，完整 strong-coupling adapter／paired commit／checkpoint 仍待完成。

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

## 下一批：移動守恆診斷（基礎驗證通過，應用驗收未完成）

Moving operator 新增 MaterialConservationDiagnostics，沿用 owned volume cells
與 required-state halo，計算各 boundary label 的 material velocity flux、wall-only
flux、原始 absolute flux sum／term count。MPI reduction 後組成與 serial
endpoint 相同定義的 wall-relative leakage、divergence theorem defect 及 discrete
moving-wall continuity，並檢查有限性與既有 roundoff identity。沒有 gather
完整數值場；目前另讀一次 halo 取得原始 flux summation 資料。

Runtime 接上 source／target geometry、publication、clock 與 audited closed-surface
volume，計算 backward-Euler volume rate、Reynolds／moving-mass defects。
PrepareCommit 在發布 prepared state 前準備診斷紀錄，FinalizeCommit 只交換
record pointers；因此 committed transition 的診斷不依賴已釋放的 source epoch。
這些是診斷與有限性檢查，尚未等同於應用層 physical acceptance gates。

測試證據目錄為 outputs/hpc03/moving-conservation-v1。第一版 endpoint 單 rank
通過；雙 rank 比較失敗於 normalized wall leakage：actual
2.9185253349559721、serial 2.9185253348966183，先前原始通量 gates 已通過。
原 logs／binary 保留；新版維持 raw-flux comparison gate，將其容許誤差連同
實際分母差傳播至 quotient comparison，避免小分母放大 roundoff 後被誤判為
原始通量錯誤。另加入 rank-local nonfinite trial 的共同拒絕與恢復後 retry。
新版 endpoint 的 1／2／4 ranks flow，以及 4 ranks pressure／traction／closed／
inertial 共 23 份 reports 全部通過，原始 physical diagnostic 最大差
1.06803e-13；source／binary／log hashes 核對於 endpoint-audit.json。單 rank
runtime 的兩步 serial conservation 對照亦通過，最大 scaled error 6.16798e-14，
場誤差 2.71512e-13，Newton／KSP 各六次；prepared 到 committed 的診斷紀錄
保持，證據為 runtime-single-audit.json。後續 2／4 ranks 也通過，最大守恆
scaled error 分別為 7.98597e-14／3.08434e-14；三種 rank counts 共七份 reports
及 source／binary／logs hashes 核對於 runtime-audit.json。這些比較驗證與串行
相同診斷定義，尚未宣稱各診斷值通過應用層 physical acceptance gates。

另加入一軸 affine contraction 的 endpoint fixture：單位截面在 0.125 秒縮短
0.03 公尺，同時保留背景分區與 active-set 改變。除串行診斷對照外，要求
closed-surface audited volume rate 及 material boundary flux 都符合獨立解析值
−0.24 m³/s。1／2／4 ranks 共七份 reports 全部通過，解析 rate 最大差
1.22125e-15，來源與 logs hashes 核對於 deforming-audit.json。此案例以任意
trial fields 驗證診斷，不屬於已收斂的變形流體守恆驗收。

## Moving graph adapter（驗證中）

ThreeDImmersedMovingDistributedFlowDomain 共用既有 graph 交易介面。
ImmersedMovingDistributedGraphBackend 持有 borrowed moving runtime 與 local-only
geometry provider；每次 coupling trial 從 accepted geometry 建立新 candidate，
共同核對 material interval 與 graph step。Provider 不能執行 MPI 或數值求解。
Graph 的 committed-field getter 始終使用 committed row ownership，不會因
trial active layout 改變而切換；提交後才讀取新 epoch 的 rows。

Port controls 先共同驗證，再寫入下一個 epoch 的 flow options；既有 immutable
committed conservation record 保留。Rollback／abort 會丟棄完整 target epoch
並還原 committed controls。Close 後 diagnostics 不再宣告 prepared，避免
late graph finalize 發布 metadata。

Backend 要求 caller 明確指定五個 finite、nonnegative normalized limits：
divergence theorem、Reynolds、moving mass、wall-relative leakage、discrete
continuity；constructor 共同檢查 policy 一致。Solve 後及 Prepare 前都檢查
當前 trial 診斷，未通過時不發布 graph measurements／committed state。
這些 limits 與 flow options 的 reference flow 必須依案例物理尺度設定，
不是跨案例通用的推薦值。

目前測試使用單位 cube、壁速 2.96 m/s、z 向 ports ±1e-4 m³/s、reference
flow 1 m³/s，兩步各平移 0.37 公尺並改變 active nodes。首次以原始 cut-cell
quadrature 執行，被 moving-mass gate 擋下（0.039495552576897464 > 0.01）；
失敗 logs／binary 保留於 outputs/hpc03/moving-graph-v1/mpi-1。改用已存在的
positive moment-fitted quadrature 後重跑，沒有放寬任何 gate。五個 limits
依序為 1e-2／1e-8／1e-2／1e-2／1e-8。

另一個收縮案例維持淨 port flow 為零，明確要求失敗源自 physical gate。
雙 rank 測試通過：moving-mass 值 0.18385222769841331 > 0.01，共同拒絕，
accepted fields／clock／graph metadata 不變，controls 還原；改回相容平移
geometry 後，同一 graph step 健康重試並成功提交。兩份 reports 與 source／
binary／logs hashes 核對於 moving-graph-v1/mass-rejection-audit.json。

完整 1／2／4 ranks 的兩步 graph／串行 field 與 port 對照仍在執行，亦包含
geometry failure、wrong material interval、force freeze failure、prepare failure、
rollback、prepared abort 與 exact accepted isolation。新 adapter 尚未接入 case
factory／正式 executor 驗收，FSI surface-field adapter 與 paired commit 亦待完成；
HPC-03D 保持未勾選。

Fitted graph 的首輪 1／2 ranks 對照在 serial oracle 被 wall fraction diagnostic
擋下，4 ranks 未啟動。獨立 probe 找到 10 個 cell 的 estimated fraction 比
certified interval 端點相差最多 9.99201e-16：Nitsche cell 層原本允許
2e-12，但 serial transient aggregation 卻要求嚴格排序。已讓 local diagnostic
記錄其既有 ordering tolerance，aggregation 累積該 tolerance 與加總 roundoff；
保留原始 bounds／estimate，沒有 clamp、重算權重或放寬 physical gates。
Expanded／compact focused regression 均通過，aggregate allowance 為
6.04601e-11，raw catalog 未改變；證據為
moving-graph-v1/fraction-wall-audit.json。原 fitted-* 失敗 binary／logs 保留。
修正後的 fraction-graph-test 已啟動 1／2／4 ranks 與 physical rejection／retry
回歸，尚未取得完整結果。

fraction-graph-test 的第一步 1／2 ranks 已完成 field／port serial 比較及所有
physical gates：最大 field scaled L2 分別為 3.11424e-13／4.23861e-13，
moving-mass 約 9.42399e-11／9.42397e-11。兩者第二步均在 geometry provider
遇到 fitted candidate cap，未完成第二次提交；4 ranks 與後續新版 rejection
重跑尚未啟動。原 fraction-* reports 保留為未通過的完整回歸。
目前先以獨立 second-endpoint geometry probe 檢查較大的明確資源上限：
max_columns=32768、NNLS max_workspace_bytes=256 MiB，candidate orders 與
moment residual tolerance 維持原值。Probe 尚在執行；未改小位移、active set
或 physical gates 以避開此案例。

較大容量 probe 已終止：繞過 candidate cap 後仍無法在既有 orders 找到正權重、
符合 moment residual 的規則，故不能只增加記憶體上限。下一個 probe 回到
原 candidate／workspace caps，將候選 support expansion 從 1.25 調為 2.0，
檢查是否為 seed support 範圍不足；所有候選仍必須位於實際流體域，moment
與 positivity gates 不變。此 probe 已完成，仍未找到符合 moment residual 的規則（細節見下）。


Support 2 probe 的失敗 cell 為 x=[1.8,2.1]、y/z=[0,0.4]；seed 的 x
範圍約 [0.01736,0.08250]，擴張後仍未覆蓋真實流體區間 [0,0.13333]。
Support 3 又揭露另一個問題：先在超出 cell 的區間放置 Gauss nodes，再丟掉
越界點，會留下不足以表示 degree-six moments 的節點間隙。現在先將候選
區間與 cell 相交，再放置 Gauss nodes；moment conditioning frame 不變。
單獨的完整第二步 geometry probe 已通過，保持原 max_columns=8192、
workspace=128 MiB、point-query cap=65536 及原 moment tolerance。

較密的 clipped grid 在非凸案例會超過候選矩陣容量；現在使用 deterministic
bounded reservoir 保留空間分布的候選子集，並恢復 parametric lexicographic
排序。每個 proposal 仍計入 point-query cap，矩陣列數仍受 max_columns 限制；
接受規則仍必須通過正權重與全部 343 個 moments 的獨立 audit。首輪 reservoir
回歸因 compact 要求嚴格點排序而失敗，失敗 log 保留；排序修正後的
ordered-fitting-regression.log 通過非凸解析 moments、原資源上限／耗盡拒絕、
expanded／compact 點值與 byte-stream 一致性。另加入實際失敗 corner cell 的
343 個解析 moments 回歸，並未使用 seed weights 當作驗收答案。

NNLS 同時修正一個獨立的 premature-stop 問題：normalized columns 的 reduced
 gradient roundoff floor 現在依 residual norm 縮放，並在每次迭代先檢查 residual
目標。解析兩欄案例 A=[(1,0),(1,δ)]、b=(1,δ/2)、δ=1e-10，舊版在
residual=5e-11 提前停止，雖然正解為 (0.5,0.5)。新版 δ=1e-8／1e-9／1e-10
及既有 scaling、random KKT、rank／capacity 回歸均通過；QR rank gate、
iteration／workspace caps 和最後的原矩陣 residual audit 不變。

ordered-graph-test 固定版本正在執行兩步 1／2／4 ranks 測試；案例明確設定
support_expansion=3，library default 仍為 1.25。新增真正的 per-field relative
L2 gate：reference norm > 1e-10 時要求 relative L2 < 1e-6；近零場則要求
absolute L2 < 1e-10。原 scaled-field 與 physical gates 也保留。此時尚未取得
完整兩步 MPI 結果，不據此勾選 HPC-03D。mesh-test 已通過。


相同 ordered binary 的 wall regression 已完成：expanded／compact 各 14 個
fraction estimate 超出 certified endpoint，最大 gap 5.55112e-16；aggregate
allowance 6.04601e-11。原始 catalog 不變，兩種模式均通過，source／binary／
rank log hashes 見 ordered-wall-audit.json。這是新版候選／NNLS 的證據，
前述 fraction-wall-audit.json 保留舊版結果。

Case factory 前置：SurfaceReaders 新增 ReadMaterialVtpBuffer／Text／Path，
使用相同有界 VTP parser 和封閉幾何驗證，但保留輸入 vertex 順序、directed
triangle connectivity 與 labels；缺省 label 解析為配置的 default_boundary_id。
拒絕 scaling／welding，避免改變材料座標或身分。ASCII、inline／appended
binary、旋轉造成 canonical 排序改變時的逐材料點插值／速度、開放表面／
無效 connectivity 拒絕，以及既有 VTP 回歸均已通過。正式 motion sidecar、
factory 所有權及 executor 接入仍待完成，不能將此 reader 當作 factory 完成。

新版 ordered-graph-test 的雙 rank physical rejection／retry 已通過：moving-mass
defect 0.1838522264663674 > 0.01，兩 rank 共同拒絕，accepted fields／clock／
graph metadata 保持不變；同一步改回平移幾何後成功提交。兩份 reports、
固定 source／binary 與 logs hashes 見 ordered-mass-rejection-audit.json。


ordered-graph-test 的完整兩步 1／2 ranks 已通過：最大 per-field relative L2
分別為 1.41086e-10／1.28498e-10，port maximum absolute error 分別為
3.28622e-16／5.59015e-16；兩次 active set 改變、兩次 accepted publication，
以及原 physical gates 均通過。4 ranks 正在執行；尚未將這一批宣告為
完整 1／2／4 ranks 驗收。

正式 ImmersedFlowCase 現已加入 prescribed_motion 分支，持有 frames、初始
MovingCutGeometry、MPI moving runtime、backend 與 graph adapter，按借用關係
順序銷毀。Provider 只做 local geometry 建立，與原 collective failure stage
配合。所有 frames 的 source coordinates／connectivity／labels／times 建立
一致性 hash，在 MPI 初始化時比對；frame file 必須是 case 內的 regular file。
明確要求 backward Euler、2–4096 個 frames、首 frame 為 time=0 的既有 surface，
frames 覆蓋設定的 simulation duration，內部 frame 邊界對齊 timestep。
每一步的 extension band 由 runtime 檢查，不能把多個 solver steps 的整個
frame interval 當作單一步位移上限。

prescribed_motion 的 schema：必填 frames（每项 time_s／surface）、
extension_layers、conservation_limits（divergence_theorem／reynolds／
moving_mass／wall_relative_leakage／discrete_continuity 五項）；可選
volume_fitting（support_expansion／candidate_orders／max_columns／
max_workspace_bytes／max_point_queries）。所有 unknown keys 拒絕，fitting
的 positivity／moment acceptance 不提供配置放寬。固定幾何的配置介面保留。

Factory 的 Distribution／SolverConfiguration 取 committed epoch；關閉前在
collective local stage 保存 geometry audit metadata，Close 後仍可讀取。
正式 CLI 的 accepted-step diagnostics 分別處理 moving epoch 的 commit_count
與全域 clock.index，不把每個 epoch 的單次 commit 誤當成整段模擬的累計值。
首輪 factory unit 已通過既有 steady／fixed transient 與 moving 載入、錯誤
配置拒絕、Begin／Abort、重複初始化拒絕及 Close 後 audit retention。它不是
完整 moving executor 求解驗收；同位移兩步的 1D–moving-3D–1D CLI fixture
已建立，正在完成新版 binary 建置後執行。

額外修正 graph interval check：驗證 material 的 start／end／evaluated time
與 graph step endpoints 相同即可，不能要求兩端相減的 DtS bitwise 等於
nominal dt。加入 start=0.1、dt=0.02 的 subtraction-roundoff fixture；錯誤
material interval 仍拒絕。此修正晚於 ordered binaries，後續以新版回歸驗證。

正式 moving CLI 首輪兩次在 preflight 拒絕，均為新 fixture 建立錯誤：
更新 time object 時漏掉子域原必填 output_every，以及只更新 3D 而未同步
1D 的 density／dynamic_viscosity。失敗 fixture／logs 分別保留；目前已回復
原 time 額外欄位，三域統一 rho=1、mu=0.1。moving-cli-2-fluid 正在使用
相同 0.37/step 位移、兩次 active set 改變的 case 執行，初始流場為正式
case 的零場，並非 graph unit 的 2.96 m/s seed。結果尚未取得。

使用方式見 [MOVING_IMMERSED_CASE.md](../MOVING_IMMERSED_CASE.md)。
下一步還須驗證一般 decimal dt 的 frame 邊界表示：例如 JSON 0.3 與
3*0.1 的一個 ULP 差異，目前仍可能被 frame coverage／alignment 嚴格比較
拒絕；這與已修正的 graph DtS subtraction check 是不同層的問題。

ordered-graph-test 的 1／2／4 ranks 全部完成，7 份 reports 的 source／binary／
logs hashes 已核對於 ordered-normal-audit.json。4 ranks 最大 field relative L2
9.14943e-11、port max abs 7.31945e-16，兩步物理門檻、rollback／fault injection
與 active-set 改變後 publication 均通過。新 factory 最終版本 unit 亦通過，
證據為 moving-final-factory-audit.json；它只驗證載入／guard／生命週期，
不代表正式 moving executor 已通過。

正式零初始場、flow inlet／pressure outlet 的 CLI 第一個 trial 在 moving-mass
gate 拒絕：0.012112992323455807 > 0.01。case、binary、雙 rank logs 保留在
moving-cli-2-fluid；未提交任何 accepted graph step。正在以另一份 fixture
啟用既有 wall_inertial_gamma0=1 的 backward-Euler 壁面阻抗，隔離檢查啟動
暫態的弱壁面條件；grid、rho／mu、位移、初始零場及所有守恆門檻不變。
新 probe 位於 moving-cli-2-inertial，尚無結論。

小數 interval 與新版 physical rejection／retry 已通過，見 duration-mass-audit.json：
start=0.1、nominal dt=0.02 的實際 subtraction roundoff 被接受，wrong interval
仍拒絕；隨後收縮案例拒絕與健康重試成功。此測試未解決 frame 邊界／
coverage 的 decimal rounding 問題，後者仍需在 factory 與 motion evaluator
共同處理，以符合正式 executor 逐次加 dt 的 clock，不能只把 frame time
改成 index*dt 就宣称解決。

Frame roundoff 修正：factory alignment／coverage 現允許至多 8*epsilon 的
相對時間誤差，並以 dt/4 限制；PrescribedSurfaceMotion 只對插值查詢做
有界 frame snap，保留 input frame times 與 published start／end／evaluated
clock 原始 bits。Snap allowance 同時受鄰接 frame duration/4 限制，不合併
過近的物理 knots；若 snapping 會把可表示的短步縮成零，回到嚴格原時間。

Prescribed motion 全部既有測試及新增十次加 0.1、0.3／0.6 速度轉折、末端
frame、實際跨 knot 拒絕、極短步及相隔一個 ULP 的 knots 均通過；既有
midpoint identity hash 未變。make mesh-test 通過。加入 dt=0.1／steps=3／
末 frame=0.3 的 factory 載入與真正不足 coverage 的拒絕回歸，正在執行。
此測試不代表已證明任意極長步數的累積 clock drift；large-step-count 驗證
仍應檢查實際累加誤差是否超過此有界 allowance。

既有 fixed-relaxation 的 inertial CLI 對照仍運作；另以完全相同物理 fixture
啟動正式 Aitken coupling 路徑，位於 moving-cli-2-aitken。這些執行使用
roundoff 修正之前的 frozen moving-final-cli（案例 dt=0.125 為精確二進位），
不把它們視作新版 decimal-boundary 的 MPI 證據。

新版 factory decimal coverage／alignment 及既有載入／拒絕／Close 回歸
均已通過，source／binary／logs hashes 核對於 roundoff-audit.json。

長 clock 累加亦加入回歸：PrescribedClockRoundoffAllowance(n) 使用
gamma_n=n*epsilon/(1-n*epsilon) 加 8*epsilon comparison margin，factory
依已知總 steps 設定 motion evaluator allowance。這是對非負 clock 逐次
加 dt 的 roundoff 上界；仍以鄰接 frame duration/4 限制 snapping，不放寬
物理門檻，也不改寫輸入 frame 或 published clock。無效／過大 allowance
明確拒絕。Factory frame 配置的 grid-alignment 比較仍只需固定的有限次
算術 margin，與執行時累積 clock allowance 分開。

一萬次加 0.1 的末時刻超過 1000：舊固定 allowance 明確拒絕；依步數的
allowance 通過，endpoint positions 與最後 frame bitwise 相同、actual clock
保留，超出 allowance 的 1000.00001 仍拒絕。完整 prescribed motion suite
通過，見 accumulated-clock-test.log；新版 factory 建置／回歸仍在進行。


固定鬆弛的 inertia CLI probe (`moving-cli-2-inertial`) 已以 exit 1 結束，
rank 0 wall time 2158.28 s，錯誤同樣為 Newton backtracking failed；這不是
通過或效能基準。新增回溯診斷保留 best candidate residual、候選數量、
最後 damping、update norm 與 linear relative residual。`backtracking-cli`
凍結版本的 Aitken probe 僅把 minimum damping 從 1/128 降到 1/8192，
物理參數與守恆門檻不變，結果另記於 `moving-cli-2-damped`，尚待驗收。


`moving-cli-2-damped` 已以 exit 1 結束，第一步接受、第二步在 Newton
iteration 27 失敗。Residual 從 59.7475554 到 23.7550261，但 14 個回溯
候選的最佳 residual 仍為 23.7583565；最後 damping=1/8192，update norm
1.730703e6，linear relative residual=3.983024e-12。減小 damping 下限
沒有完成原案例，詳見 `backtracking-failure-audit.json`。後續使用同一
`backtracking-cli`、同物理條件／幾何路徑／守恆門檻，僅將 dt=.125、2 步
改為 dt=.0625、4 步（結束時間仍 .25），診斷時間解析度影響。
新 fixture 與輸入 hashes 為 `cli-fixture-half-step*`，作業為
`moving-cli-2-half-step`；即使新案例通過，也不代表原 dt=.125 案例通過。


`moving-cli-2-half-step` 在第一個 target geometry 建構即以 exit 1 結束：
`fitted cut-cell support has zero width`，未執行 Newton。兩個 rank reports
及錯誤存於 `half-step-failure-audit.json`。目前 fitting 使用 seed point
各軸 min/max 建立 conditioning frame；只要某軸 seed 座標相同便拒絕。
後續需確認薄切割體積的 seed 是否合法但不足以估計 support，並補上適當
候選範圍與完整 degree-six moment 驗證；此診斷不能支持減小 dt 後的
Newton 成功或失敗結論。


薄切割 seed 的直接探查確認：dt=.0625 第一個 target 的 cell 0 有 144
個 seed points，x 坐標全部為 .9826420389492565。對應真實 box intersection
為局部 [.95,1]×[.25,1]×[.25,1]，正體積，故 seed 零跨度不是空體積證明。
最初只改用 full-cell conditioning 雖通過單點 seed 的非凸 343 moments
測試，實際薄切割仍 NNLS exhausted（relative residual 約 1.277e-4），
原 binary／header／失敗 log 已保留為 `full-cell-spread-*` 及
`zero-spread-fitted-probe.log`。

現改為只在 seed 零跨度時，裁切 closed surface triangles 至 cell，並加入
Inside／Boundary 的 box corners，重建未解析方向的支撐範圍。支撐 extrema
包括裁切邊界交點與保留的 box corners；非退化 seed 沿用原路徑。新增的
八次 corner classification 計入原 point-query cap，triangle 數沿用 moment
上限；ambiguous corner 仍拒絕。候選 selection、NNLS／完整 343 moments
門檻與正權重要求不變。新增局部 [.95,1] 薄切割的解析 moments 回歸，
並重新驗證完整 target geometry；結果另記 `clipped-spread-*`。


裁切支撐範圍修正的 unit 與完整第一個半步 target geometry probe 均已
exit 0。單點 seed 非凸體積、x 零跨度薄切割體積各 343 個解析 moments、
正權重及 query cap 拒絕，連同原 compact／expanded／diagnostics 回歸通過；
完整目標幾何在原 fitting caps 下成功。來源／binary／logs 雜湊與驗收
範圍見 `clipped-spread-audit.json`。Probe 來源的 evaluated time 是 .0625
（舊輸出 label 仍寫 second_endpoint，不代表 .25 的端點）；完整 CLI
重建後另行執行，幾何 probe 不取代 Newton 或 graph 驗收。


修正後完整雙 rank CLI (`moving-cli-2-half-step-spread`) 已 exit 0、無 timeout，
接受四步至 t=.25。Source／fixture／binary／兩 rank logs／全部 CSV hashes
已核對於 `half-step-spread-cli-2-audit.json`。Graph 迭代數為 7／4／4／4，
八個 accepted edge records 的 normalized pressure／flow residual 均符合
1e-6／1e-10。最大 moving mass 與 wall-relative leakage 約 .002728381，
最大 Reynolds defect 3.72e-15；四個 accepted Newton trial 的 minimum
accepted damping 均為 1。原 dt=.125 的第二步失敗仍保留，不以此覆蓋。

Rank 0 wall time 1283.76 s；profile exclusive geometry 820.38 s、assembly
433.43 s、solver setup 9.60 s、linear solve .156 s，所有 rank RSS 保留在
run reports。這些是在其他回歸同時執行時取得，不是獨立效能基準。
同一凍結 CLI／fixture 的 1／4-rank sweep 已啟動，另存
`moving-cli-{1,4}-half-step-spread`，後續比较完整 port／graph history。

新增 `scripts/hpc_make_moving_graph_case.py`，由既有 repository chain fixture
建立 prescribed moving case 與 SHA-256 manifest，不自動求解、不覆寫既有
目錄。預設四步與 `--steps 2` 的全部輸入檔案，分別與已驗收 half-step
fixture／已記錄失敗的 coarse-step fixture 逐 byte 相同；既有目錄拒絕後
所有檔案 hashes 不變。驗證於 `moving-fixture-reproduction-audit.json`，
使用命令已加入 `docs/MOVING_IMMERSED_CASE.md`，因此不必依賴 outputs
中手工修改的案例才能重現這批驗收。


保存結果的比較另交付 `scripts/hpc_check_moving_graph.py`：逐 rank 檢查
退出／timeout／log digest、row coverage、accepted clock、CLI 已記錄的三個
moving physical diagnostics、accepted graph edge gates 與完整 port coverage；
可指定 reference run 沿用既有 graph port tolerance。JSON 明確列出限制，
不以此取代 full fields／source provenance／未輸出的 physical gates 或效能
驗收。正常雙 rank 結果與七種錯誤注入已通過（缺少／失敗 rank、損壞 log、
錯誤 clock、缺 port、NaN、錯誤 reference），證據於
`moving-checker-negative-audit.json`。重現命令改用 rank wrapper 保存完整報告。

同一凍結 CLI 的單 rank 四步作業現已 exit 0、未 timeout。保存結果 checker
以單 rank 為 reference 比較 1／2 ranks，全部 rank／clock／graph gates／port
coverage 檢查通過；雙 rank 最大 port error／tolerance 比值為
3.239808816846962e-05（門檻為 1），沿用 `1e-12 + 1e-6 * abs(reference)`。
結果及 logs／CSV hashes 保存於 `half-step-cli-1-2-comparison.json`；
4 ranks 作業仍在執行。這是相同輸入的端口歷程比較，不是完整場比較或
獨占 scaling 證據。

正式 CLI 的 4 ranks 作業亦已完成，七份 rank reports 全部 exit 0、未
timeout。`half-step-cli-1-2-4-comparison.json` 驗證三種 rank 數的四步
accepted clocks、graph gates、port coverage 與對單 rank 的端口數值門檻。
這完成此 prescribed-motion 小案例的本機 MPI matrix。

後續 moving-fluid adapter 已接入分散式 surface capture、material composition、
traction publication、全域守恆 gate、single-owner bounded membrane、paired commit
與 checkpoint。非零強耦合通過 1／2／4 ranks；moving flow 4→2／4→1 與完整
FSI pair 4→2、4→1、1→4-rank 新作業續算也通過 stable-ID ownership 搬移。
因此 active-set 更新、history extension、提交／回復、正式 moving CLI 與
restart 已滿足 HPC-03D。現有 profiling 顯示幾何建構主導，但沒有支持執行中
重新分配固定背景的穩定效益，故保留固定背景 partition。大型記憶體、正式
moving／FSI native graph CLI 與跨節點要求由 HPC-05C、HPC-06、HPC-07D／09
追蹤；不影響 HPC-03D 的完成狀態。重分區證據見
[HPC-05D 完成報告](HPC_05D_MOVING_FSI_RESTART_REPORT.md)。
