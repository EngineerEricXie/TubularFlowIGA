# HPC-07A：owned 表面 publication 到封閉材料幾何

狀態：材料合成元件、moving runtime 的 owned coefficient capture 與 stamped
traction publication 已通過 1／2／4 ranks。Same-trial capture 與 accepted-state
隔離保持 exact；重新求解以既有逐場 scaled L2 門檻驗證。延後資源釋放的
local abort 入口另已通過 4 ranks 驗收。Moving-fluid FSI adapter 已接上真實
膜 runtime 與 paired commit：全 clamped 交易通過 1／2／4 ranks，非零
變形交易與十進位跨步各通過 4 ranks。非零案例的 owned 場值矩陣已通過
1／2／4 ranks；非零強耦合兩步及完整 field／history／port／守恆比較也已
通過 1／2／4 ranks，提交前故障與迭代耗盡 rollback／retry 已分別驗證。
膜 checkpoint 支援跨 rank 數還原；moving fluid 的材料／owned field／
守恆紀錄已接到實際 MPI 檔案 bundle，4-rank changing-layout 續算通過。
完整 FSI paired restart 已通過同 rank、新作業及 4→2、4→1、1→4-rank
重分區續算；HPC-03D、HPC-05D、HPC-07A／B 已勾選完成。全域強耦合與
完整 restart 的跨節點驗收仍由 HPC-07C／D、HPC-09 追蹤。

每批凍結 source／binary 的精確對應另核對於
`outputs/hpc03/moving-graph-v1/batch-source-resolution-audit.json`。

目前 paired runtime 驗收可由 repository 直接重現：

```bash
make -C solvers/cpu immersed_moving_distributed_fsi_runtime_test \
  CXX=mpicxx PETSC_DIR=/path/to/petsc PETSC_ARCH=your-arch
mkdir -p outputs/fsi-validation/mpi-4
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -n 4 \
  python3 scripts/hpc_rank_run.py --output-dir outputs/fsi-validation/mpi-4 \
  --expected-ranks 4 --timeout 5400 -- \
  solvers/cpu/immersed_moving_distributed_fsi_runtime_test deforming \
  -immersed_moving_mat_mumps_icntl_14 100
```

省略 `deforming` 可執行全 clamped rollback 與 decimal-step 案例。非零
案例須以同一 binary 分別執行 1／2／4 ranks，使用各自全新的輸出目錄，
再比較保存結果：

```bash
python3 scripts/hpc_check_moving_fsi.py \
  --reference-run outputs/fsi-validation/mpi-1 \
  --run outputs/fsi-validation/mpi-1 \
  --run outputs/fsi-validation/mpi-2 \
  --run outputs/fsi-validation/mpi-4
```

Checker 檢查每份退出報告與 log digest、完整唯一 owned rows／surface IDs，
並分量比較速度、壓力、輔助未知數與十二個 surface fields。正常合成輸出
及 11 種損壞注入已通過，見 `paired-fsi-checker-negative-audit.json`；
這項工具測試不代替實際求解矩陣，也不驗證 source provenance 或效能。
以下歷次結果各有自己的驗收範圍，不能合併解讀成同一 HEAD 的完整驗收。

`MaterialSurfacePatchKinematics::ComposeCompletePatch` 是共用的 geometry-only
核心：要求完整、依 map 排序的所有節點，檢查有限值、夾持 seam exact zero、
committed material／topology／clock、backward-Euler 位移速度一致性與封閉表面。
既有 serial `ComposeTarget` 仍先拒絕非完整單分區 publication；不讓 partial
slice 經由此 API 靜默變成完整幾何。舊／新 header 各自編譯的獨立 probe
確認七個 target／material／topology／composition identities 完全一致。

`ComposeDistributedMaterialSurfacePatch` 先驗證外部 transaction authority
提供的 expected stamp／committed identity、共同 reference／context 與資源上限，
再經 existing owner-to-ghost 通道取得閉合幾何需要的 patch coordinates。
來源 publication 保持 owned，函式不提交狀態，也不複製或求解完整 fluid field。
輸出 composition identity 綁定所有 rank 的來源 publication identities；最後
要求目標 content identity 在 communicator 中一致。

目前封閉幾何仍在每個 rank 複製，因此這個橋接需要把完整 patch kinematics
送到每個幾何 owner。預設上限為 4096 patch nodes、1000000 full material
vertices，還受 PointIdentityLimits 約束。這不是 distributed geometry storage，
也不是非匹配介面的投影實作。

回歸：兩步非零 patch displacement、對串行 target identity 比較、原 input
保持不變；2／4 ranks 中 rank 0 為空 owner。錯誤 expected stamp、stale
committed identity、夾持 seam 非零、BE 不一致、NaN、缺少 node owner、
資源上限拒絕均共同完成；之後健康重試產生相同 composition identity。
7 份 reports 的 source／binary／logs hashes 核對於
`outputs/hpc07/material-composition-v1/mpi-audit.json`；舊新 identities 證據在
`legacy-identity-audit.json`。既有 serial material patch suite 也通過。

```bash
make -C solvers/cpu distributed_material_surface_patch_test \
  material_surface_patch_kinematics_test PETSC_DIR=/path/to/petsc
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  mpiexec -n 4 solvers/cpu/distributed_material_surface_patch_test
```


移動幾何的 traction provenance 修正：canonical triangle 排序會隨座標改變，
因此不能用 reference geometry 的 canonical ID 選取 current quadrature。
`MaterialSurfacePatchMap::LayoutTrianglesByCanonical(current)` 現在依穩定的
source triangle ID 建立當前映射，先驗證 material／topology identity；串行
traction、owned traction extraction 與既有 moving FSI 的 surface-state
選取都使用這個映射。原 reference lookup API 與檔案格式不變。

旋轉封閉 tetrahedron、保留原 patch map 的回歸確認 canonical 排序確實改變，
串行及 owned extraction 的常壓 traction 與合力仍符合解析值。material map
suite 與既有 traction suite 的 1／2／4 ranks 均通過，紀錄於
`outputs/hpc03/moving-graph-v1/current-provenance-traction-audit.json`。
旋轉載荷檢查目前是每個 rank 的 local extraction；分散式載荷投影沿用
原始幾何測試，不能據此宣稱完整 distributed moving FSI 已驗收。


下一段接線已加入 moving operator 的 `CaptureOwnedPatchState` 與 runtime 的
`CaptureTrialPatchState`：按當前 material triangle provenance 選取 owned cut
cells，經現有 required-state halo 讀取四個流體係數，不 gather 完整流場。
Runtime 僅允許 converged、尚未 prepare 的 trial 擷取，避免發布初始 seed
或已準備提交後的可變資料。這仍是 traction producer 的輸入橋接，尚未
接通完整 fluid／membrane transaction。

擷取測試正在補齊：operator 對既有完整小案例係數逐值比較 owned cell
coverage 與所有 nodal coefficients；runtime 測試未開始／未收斂／prepared／
committed 階段拒絕，以及 abort／retry 擷取一致性。首輪測試因將帶非零
運動速度的 initial evaluation 誤作靜止 reference 而拒絕，已保留
`patch-capture-4.log`；修正測試為相同材料 reference vertices、零速度後，
重建並另存 `patch-capture-reference-*`，不放寬 reference-state guard。


既有完整 moving FSI runtime 回歸在 600 s timeout 結束（exit 124、未產生
完成數值報告），不是數值通過證據。相同凍結 binary 已改用 3600 s 時限
重跑，紀錄於 `current-provenance-fsi-retry.log`；原逾時 metadata 與 hashes
保留在 `current-provenance-fsi-timeout.json`。Operator 擷取的 1／4 ranks
已通過精確係數比較與原殘差／Jacobian action／守恆門檻，2 ranks 與
正式 runtime 擷取生命週期測試尚待完成。


Operator 擷取的 1／2／4 ranks 已全部通過，七行 rank results 及 hashes
核對於 `patch-capture-operator-audit.json`；最大 residual relative L2 為
3.11e-16、Jacobian action relative L2 為 2.20e-15，最大原守恆診斷比較
誤差為 1.069e-13。正式 runtime 首輪已通過新增的 solve／abort／retry
擷取段落，但隨後 serial oracle 因未提供既有 MUMPS 設定而 KSP 失敗，
整體 exit 1，不能記為通過。已按 `moving-conservation-v1/runtime-4` 的
原驗收 argv 補回 oracle MUMPS／GMRES 與 workspace 設定，重跑兩次
active-layout 改變的 `changing` 案例，紀錄在
`patch-capture-runtime-changing-2.log`；物理門檻未變。


正式 moving runtime 已加入 `BuildTrialSurfaceTraction` 接線。它先核對
caller 的 FSI context（step／start／dt／end）與 expected material content，
再擷取 converged、unprepared trial 的 owned patch coefficients，從該
runtime capture 建立 producer identity。使用與 volume assembly 相同的
固定背景分區（包括 inactive cells）送入既有分散式 traction publication，
保留合力／力矩／離散功 gate 與 bounded single-owner projection。
這個入口不提交 fluid 或 membrane state；完整強耦合 transaction 仍待接線。

新增驗收預先固定為兩次 active-layout 改變後，owned nodal force 與
traction 對獨立串行 flow oracle 各分量絕對誤差 <1e-8；此平移封閉流體
案例的解析 stress 接近零，非零壓力／黏性解析 traction 仍由既有
traction suite 覆蓋。另測單一 rank 錯誤 step／material identity 的共同
拒絕與健康發布。這些新 assertions 的執行結果尚待取得。


`patch-capture-runtime-changing-2` 現已 exit 0：兩步 active-layout 改變、
未開始／未收斂／prepared／committed 擷取拒絕、converged 擷取與
abort／retry 係數逐位元相同均通過。原流場對 serial oracle 的最大 scaled
L2=4.05808e-13、守恆比較最大 scaled error=7.98597e-14，Newton／KSP
合計各 6 次。來源、凍結 binary 與 log 的精確 hashes 已核對於
`patch-capture-runtime-audit.json`。此 binary 早於新的 traction publication
入口；該入口另外由 `traction-publication-runtime-test` 凍結版本、
`traction-publication-2` 作業驗收，不混用兩批證據。


`traction-publication-2` 已通過（兩份 rank reports exit 0、無 timeout）：兩步
active-layout 改變後的 owned force／traction 均符合預先固定的 serial
oracle 比較門檻；單 rank 錯誤 step／material 的共同拒絕與健康發布通過，
既有全域合力／力矩／離散功 gate 保持。來源、binary、所有 logs 與數值
assertions 的範圍核對於 `traction-publication-audit.json`。Runtime
publication 的 1／4 ranks、實際膜 feedback、共同 commit 與 checkpoint
仍未由這批測試驗收。


Runtime traction publication 的單 rank 作業已 exit 0（411.64 s），原四 rank
作業則在 1200 s timeout，四份 reports 都標記 124；沒有完成輸出，不能
宣稱通過。測試新增第一個 solve／oracle／traction publication 等階段標記，
並把 retry coefficient／publication stamp assertions 放進 collective error
boundary，避免局部測試例外使其他 ranks 等待。未改數值門檻；新的
`traction-stage-4` 使用獨立 binary／source manifest 與 1800 s timeout。
原完整 serial moving FSI 回歸已輸出 channel converged／守恆結果且程序
已結束，但工具 session 遺失，沒有持久化退出碼；此 log 只能支持數值
輸出觀察，完整退出驗收仍待正式 wrapper 重跑。


四 rank 階段診斷已有限時間 exit 1，四份 stderr 都共同定位於
`moving patch retry comparison: moving patch retry coefficients differ`，
尚未進入第一步 oracle／traction publication。兩份 run reports 來不及
在 launcher 的 nonzero abort 前發布，已於 `traction-stage-failure-audit.json`
標記缺失，沒有補造結果。下一個診斷保留 exact comparison，增加差值與
coefficient magnitude，並設定 launcher 不提前中止非零退出的 wrappers，
讓各 rank 的結束報告有機會完整發布。


四 rank 差值診斷取得四份有限時間 exit 1 reports：重新 nonlinear solve
後 max absolute coefficient difference=6.248335182590381e-13，係數最大
量級=3.8210132778151995；既有 accepted-state bitwise isolation 已通過。
`traction-retry-delta-audit.json` 保存各 rank reports／logs 核對。停用
launcher 提前 abort 後 launcher 本身回傳 0，因此後續 driver 必須逐一
檢查 rank return codes，不能只看 mpiexec exit status。

原新增測試把重新數值求解誤當作同一狀態的複製而要求 bitwise equality。
現在分開驗證：同一 converged trial 連續擷取仍 exact；accepted state／
geometry 的隔離與 publication 仍維持原 exact checks；重新 nonlinear
solve 的 patch coefficients 使用本測試既有 oracle 的逐場
`sqrt(sum(error²))/max(1,sqrt(sum(reference²))) < 1e-8`。
這項修正不更動 production solver、Newton stopping 或物理／牽引力 gate。
四 rank 的原 exact failure 保留，新的 4／1／2-rank suite 另存
`traction-retry-contract-*`，輸出實際 patch retry scaled L2，等待驗收。

另外，完整 moving FSI test 的 Makefile 依賴補上 header 集合，使間接使用的
material map／quadrature 等 header 修改也會觸發重建，避免誤跑舊 binary。


`traction-retry-contract-4` 現已四份 rank reports exit 0、無 timeout；兩步
active-layout 改變的 field／traction oracle、全域守恆、錯誤 context 拒絕
與健康發布、commit／rollback 均通過。同一 trial 重複擷取 exact；重新
求解 patch scaled L2=9.88779e-14，最大 field scaled L2=2.16864e-13，
最大守恆 scaled error=3.58394e-14。完整來源、binary 與 logs 核對於
`traction-retry-contract-4-audit.json`。原 bitwise retry assertion 失敗及
timeout 保留，不能將其改寫為成功；同一新 test revision 的 1／2 ranks
仍由 driver 接續執行。此結果尚不涵蓋完整 membrane feedback 或 checkpoint。

成對 FSI commit coordinator 的 failure gate 要求 local／noexcept rollback，
原 moving runtime 的 `AbortTrial()` 會立即 collective Close PETSc 資源，
不能直接接入。新增 `AbortTrialDeferredNoexcept()` 先撤回 trial clock 與
prepared conservation，將整個 candidate epoch 移至原有 retired slot；
BeginTrial 或 Close 再共同釋放。已接受的 geometry／Vec／clock 保持不變，
沒有新增完整場複本，最多仍保留一個 retired epoch。這個入口須在已共同
決定失敗後由所有 ranks 呼叫；本身不提供 MPI process-loss recovery。

測試新增 prepared trial 的 local abort、未 prepared failure 的 repeated
local abort、舊 accepted state 精確隔離、撤回後拒絕 trial diagnostics、
下一次 BeginTrial 清理並健康重試。新的 source revision 正在獨立建置／
驗收，不能用先前 traction-publication suite 的通過紀錄替代這個入口。


同一 `traction-retry-contract-test` 的 1／2／4 ranks 現已全數通過，七份
rank reports 與全部 log hashes 核對於
`traction-retry-contract-matrix-audit.json`。這批只涵蓋既有 transaction 與
traction 接線；新 deferred-abort 入口使用獨立
`deferred-abort-source-manifest.json`／`deferred-abort-4` 作業驗收。

`deferred-abort-4` 現已完成：四份 rank reports 全部 exit 0、未 timeout，
source／binary 與 logs hashes 核對於 `deferred-abort-audit.json`。Prepared
與未 prepared trial 的撤回、重複撤回、accepted state 精確隔離及下一次
BeginTrial 健康重試均通過；兩次 accepted steps 的最大場 scaled L2 為
3.94291e-13，最大守恆 scaled error 為 4.35693e-14。這是單工作站
runtime 回歸，不代表完整 FSI coordinator 已接線，也不是獨占效能測試。

新增 `ImmersedMovingDistributedFsiRuntime.hpp`，將 owned kinematics 的權威
stamp／committed material 檢查、完整材料合成、local geometry provider、
分散式 moving solve、五項必填物理門檻與 owned traction 發布串接。
Adapter 保存確切 FSI context（含 coupling iteration），並提供 paired
coordinator 的 local finalize／deferred abort hooks；發布與 accepted fluid
epoch 在同一 finalize 階段交換。借用的 runtime 與 patch map 必須由 adapter
獨占交易期間的修改權。這是新實作，僅完成真實膜 runtime／coordinator
template 的編譯檢查，數值與故障注入測試尚未完成，不據此宣稱完整 FSI 驗收。

另以正式 rank wrapper 完成既有串行 FSI 回歸 `current-fsi-batch-1`：exit 0、
未 timeout，9-node／8-triangle channel 收斂，open flow=-0.00054、moving
mass diagnostic=0.0244871、resultant=0.0107223。精確來源、binary、報告與
logs hashes 保存於 `current-fsi-batch-audit.json`。沿用此串行案例原本的
物理門檻，不能套用到 prescribed graph 或新分散式 adapter 的驗收。

新 adapter 的 paired regression 使用真實 moving fluid 與 single-owner
membrane runtime，注入 structure prepare 失敗及單 rank context mismatch，
檢查 accepted field／material／clock／雙方 publications 不變及健康重試。
此案例四個 patch nodes 全部 clamped，用於隔離交易語義；可變形 membrane
feedback 與強耦合收斂仍須另外驗收。首輪 `paired-fsi-4` 因測試 fixture
fluid interface 的 subsystem 名稱不一致，四 ranks 均在 preflight exit 1；
紀錄保存於 `paired-fsi-fixture-failure-audit.json`。已修正 fixture，後續
驗收使用獨立 binary／輸出，不覆寫失敗證據。

修正後 `paired-fsi-subsystem-4` 四份 rank reports 全部 exit 0、未 timeout，
source／binary 與 logs 核對於 `paired-fsi-subsystem-4-audit.json`。Rank 0
為空 surface owner；structure prepare 故障、單 rank coupling-iteration
錯誤、accepted state 精確回復、健康重試與雙方 publication 成對提交均通過。
同一凍結 binary 的 1／2 ranks 正在接續驗收。這仍是全 clamped 交易案例，
不能取代可變形表面、完整強耦合或 checkpoint 驗收。

同一 `paired-fsi-subsystem-test` 的 1／2／4 ranks 已全數通過，七份 rank
reports 與 logs hashes 核對於 `paired-fsi-subsystem-matrix-audit.json`。
下一個 `deforming` 變體沿用 `CompliantChannelFsiFixture.hpp` 的 9-node
patch、八個 clamps、一個自由中心節點、3³ 背景網格與流體／膜材料參數；
僅將初始 material clock 平移到 0，以符合目前 single-owner membrane 的
初始時鐘。首個 trial 施加中心位移／速度 .006，要求流體牽引力與膜回應
非零，再執行 prepare／context failure 與 exact rollback。事先設定 gate
為 divergence／moving mass／wall leakage 各 .03、Reynolds／discrete
continuity 各 1e-8。此變體尚待數值驗收；成對提交測試不宣稱 FSI fixed
point 已收斂，後者仍須完整 feedback 迴圈及串行／多 rank 對照。

跨步複查另發現 traction transaction guard 對 nominal `context.dt_s` 與
material endpoint subtraction 作 bitwise 比較；例如 .325 + .1 的端點
差為 .10000000000000003，會拒絕合法 context。改為沿用 graph 的端點
契約：ValidateFsiTrialContext 後仍精確核對 step、start、end 及 material
identity，不額外要求浮點相減等於名義步長。Paired regression 加入首步
.125 後連續三步 .1，逐步檢查 fluid／structure publication 與 accepted
clock 一致；新 revision 正在獨立建置驗收。已啟動的 deforming binary
仍保留原來源，不受此次修改影響。

`paired-fsi-deforming-4` 已通過四 ranks 的非零 patch force／膜位移、
prepare／context 故障、精確回復與健康提交，報告與 log hashes 見
`paired-fsi-deforming-4-audit.json`。`paired-fsi-decimal-4` 亦通過首步
.125 加三步 .1 的雙方 publication／clock 檢查，見
`paired-fsi-decimal-4-audit.json`。兩者使用各自凍結來源；前者沒有跨 rank
場值比較，後者全 clamped，均不能證明多步強耦合。

為補足數值矩陣，新 test revision 輸出 test-only owned fluid rows 與表面
位移／速度／牽引力／consistent force，不增加 production gather。
`scripts/hpc_check_moving_fsi.py` 檢查完成報告、hashes、唯一 owned coverage
及 reference IDs，沿用每個分量 `L2(error)/max(1,L2(reference)) < 1e-8`。
新 1／2／4-rank matrix 與 checker 的實際結果仍待取得。

新增 `DistributedStrongFsiStep.hpp` 將真實 adapter／membrane 的迭代串接：
以 immutable reference facets 建立 owned normals，重用全域 area-weighted
RMS 收斂與 distributed weighted Aitken；每次未收斂時雙方撤回 trial，從
相同 accepted epoch 重新求解。Relaxed 壁面速度以 accepted **fluid**
geometry 的位移差除以 dt，避免上一時間步有限容差留下的 fluid／structure
位置差使材料 BE 檢查失效。除既有位移 RMS 門檻外，也要求最大 Cartesian
速度殘差乘 dt 不超過相同位移門檻；並未放寬材料一致性檢查。

所有可失敗的 result allocation／驗證都在 paired finalize 前完成，例外或
最大迭代數耗盡時撤回雙方 trial。`strong` 測試模式沿用 channel 的材料、
physical gates 與 16-iteration Aitken 設定，執行兩個連續時間步；新迴圈
已完成 template 編譯檢查，數值、失敗注入與多 rank 比較尚待驗收。

`strong-zero` 模式另使用全 clamped 的真實 fluid／membrane runtime 隔離
新 coordinator 的例外生命週期：注入收斂後的 precommit failure、單 rank
structure prepare failure 及單 rank 無效 maximum_iterations，要求 accepted
owned field／material／clock 不變、雙方 trial 與未建立的 committed
publications 仍不可讀，再健康推進兩步。`strong-fsi-failure-test` 已獨立
編譯並啟動 `strong-fsi-failure-4`，結果尚待取得；這不取代非零 fixed-point
收斂或迭代上限耗盡的測試。

`strong-fsi-failure-4` 現已完成，四份 rank reports 全部 exit 0、未 timeout；
三種拒絕／精確回復與後續兩步健康提交均通過。兩步各一次迭代，位移 RMS
及速度殘差為零。來源、binary、report 與 log hashes 核對於
`strong-fsi-failure-4-audit.json`；非零兩步作業 `strong-fsi-4` 仍在執行。

非零 `strong-fsi-4` 已結束：四 ranks 在首輪 divergence conservation gate 0
一致 exit 1，尚無 accepted step，見 `strong-fsi-gate-failure-audit.json`。
Fixture 沿用的 flow-controller reference 為 1e-12 m³/s；moving normalization
使用此 floor、體積變化與淨表面流量的最大值，靜止壁面 through-flow 的
入出流量相消可能使尺度過小。目前僅是待核對原因，不能把失敗當作通過。
Adapter 錯誤訊息已加入實際 normalized value、limit、normalization scale
與原始 divergence／moving-mass／wall-leakage 缺陷；先以同一案例與原門檻
取得診斷，再決定應修正案例物理尺度或數值離散。

同門檻診斷重現於 `strong-fsi-gate-detail-audit.json`：normalization=
1.7289240792883777e-6 m³/s、divergence defect=-1.4822239579289541e-7 m³/s，
normalized divergence=.085731003210912246；moving mass defect 則恰好等於
normalization，因此 normalized mass=1。零初始壁速的 through-flow predictor
沒有體積變化，default 1e-12 floor 太小，實際以淨流量缺陷正規化自身。

強耦合測試現明確設定固定參考流量 1.8e-3 m³/s，由輸入資料計算
`width * height^3 * pressure_drop / (12 * viscosity * length)`，其中 width、
height、length 均 .6 m、pressure drop=.1 Pa、viscosity=1 Pa·s。這是忽略
側壁的平行板流量上界，作為此小方管的特徵流量，不宣稱是方管精確解。
其推導為積分 `u(z)=Δp*(h²/4-z²)/(2*mu*L)`。尺度不依量測殘差或 rank
數決定；比例門檻仍 .03／1e-8／.03／.03／1e-8。這會改變對應的絕對
容許缺陷，舊失敗保留；新的强耦合結果獨立驗收，未更改既有 paired
deformation matrix、流體方程、壓差或材料參數。

Owned field matrix 的單 rank 已完成：864 個 fluid unknowns 與 9 個 surface
nodes 具完整唯一覆蓋，free center displacement=-1.5172251536063165e-6 m。
`paired-fsi-fields-1-audit.json` 核對精確 source／binary 與保存輸出；目前
只完成自身格式檢查，2／4 ranks 比較仍待完成。

`strong-exhaustion` 模式新增真實非零首輪的 iteration-limit 回復：先以
maximum_iterations=1 執行，必須精確收到 fixed-point 未收斂錯誤，不能
用 geometry／solver／physical gate 拒絕冒充耗盡。拒絕後檢查 owned field、
material、clock 精確不變與雙方 publication 不可讀，再以正常 16-iteration
控制重試兩步。此變體仍沿用固定物理 reference flow；尚待執行結果。

`paired-fsi-fields-1-2-comparison.json` 已核對 1／2 ranks 的完整 owned
fluid／surface output。最大分量 scaled L2=1.5004551024690282e-12，來自
壓力；表面十二個分量的最大值為 1.2435806593839359e-14。皆低於既定
1e-8 門檻，4 ranks 仍待完成。此對照使用同一凍結 source／binary，
不是獨立 serial discretization oracle，也不是強耦合 fixed-point 結果。

Strong step 另加入進入 numerical collective 前的 communicator preflight：
流體、膜與 step communicator 必須 IDENT 或 CONGRUENT，surface partition
rank/count 必須對應；不能等到 paired finalize 才發現錯用 communicator。
`SingleOwnerMembraneRuntime` 提供只讀 communicator getter。`strong-zero`
新增在 world runtimes 上錯用 MPI_COMM_SELF 的拒絕與不變狀態檢查，
新 revision 尚待獨立驗收。

`strong-fsi-comm-4` 已通過四 ranks 的 communicator mismatch 拒絕、
precommit／structure prepare／invalid controls 故障、exact rollback 與
兩步健康重試，見 `strong-fsi-comm-4-audit.json`。九個近期 FSI 凍結版本
的全部 source／binary hashes 已重核於 `fsi-batch-source-resolution-audit.json`；
這只證明來源可追溯，不把仍執行中的非零強耦合作業視為通過。

非零成對交易的場值矩陣已全數通過：`paired-fsi-fields-matrix-comparison.json`
比較 1／2／4 ranks 的 864 fluid unknowns 與 9 個表面節點，七份 rank
reports 均 exit 0、未 timeout。最大分量 scaled L2 分別為 0、
1.5004551024690282e-12、1.340721349608626e-12；表面最大分量誤差分別為
0、1.2435806593839359e-14、1.227811051742249e-14，均低於 1e-8。
`paired-fsi-fields-matrix-audit.json` 綁定來源、binary 與比較結果 hashes。

固定物理 reference 的 `strong-fsi-flow-scale-4` 已觀察到第一步接受：
6 次迭代、位移 RMS=1.9011931674277533e-8 m、最大速度殘差乘 dt=
3.8023863348555067e-8 m、threshold=6.5417913784292131e-8 m；第二步
仍在執行，不能視為完整 job 通過。`strong-fsi-exhaustion-1` 亦已輸出
iteration exhaustion 的 exact rollback 通過，正常兩步重試尚未完成。

Checkpoint 的材料幾何元件 `MaterialSurfaceCheckpoint.hpp` 已加入：使用
既有 bounded little-endian metadata codec，保存 reference／current vertices、
velocity、穩定來源 connectivity／labels 與 evaluation／step clocks。讀回時
必須提供 bundle authority 預期的 content identity，重新建立封閉幾何並比對
完整內容身分；不僅相信 payload 內自述的 hash。Geometry 預設上限為
100000 vertices／200000 triangles，另受 metadata codec 的 16 MiB 上限
約束；這是 replicated geometry metadata，不是 fluid field shard。

`make -C solvers/cpu material_surface_checkpoint_test CXX=mpicxx PETSC_DIR=...`
與 executable 執行已通過，參考／非零變形幾何皆逐 byte 往返一致，包含
source 與 canonical provenance；截斷、trailing bytes、錯誤 authority、
labels 損壞、非有限座標、上限與 forged counts 均拒絕。來源與結果綁定於
`material-checkpoint-source-manifest.json`。尚未接上 fluid shards、膜狀態、
atomic bundle、restore candidate 或完整續跑，因此 HPC-05D／HPC-07D
仍未完成；超出 metadata 上限的大型幾何仍需分塊路徑。

`strong-fsi-flow-scale-4` 現已完整完成：四份 rank reports 全部 exit 0、
未 timeout。第一／第二步分別 6／7 次迭代；第二步位移 RMS=
3.6118189853734861e-9 m、最大速度殘差乘 dt=4.5247501270581233e-8 m，
均低於 threshold=6.6069268614094625e-8 m。精確 source／binary、七項
每步欄位與 logs hashes 核對於 `strong-fsi-flow-scale-4-audit.json`。
這完成單工作站 4-rank 非零兩步強耦合；尚未完成其他 rank 數的收斂歷程／
fields 比較、正在執行的 exhaustion retry，以及完整 checkpoint。

下一批 `strong` test revision 在每個 accepted step 輸出完整 owned field
與穩定 node ID／component mapping、表面十二分量、端口 area／flow／pressure／
normal traction、五項守恆值及其門檻，另保存每次 Aitken 的 RMS／max／
scale／threshold／velocity residual／relaxation。這些是 regression 輸出，
用於後續 1／2／4 ranks 的逐步場值與歷程對照；尚未據此宣稱矩陣通過。

`scripts/hpc_check_strong_fsi.py` 讀取上述兩步輸出，要求每個 rank 正常
退出、log hashes 一致、owned rows／surface IDs 唯一完整、node/component
映射一致、兩個 ports 與守恆紀錄完整，以及各 rank 的 shared history 完全
一致。逐步場值沿用分量 scaled L2 <1e-8；迭代、端口與守恆數值沿用
`1e-12 + 1e-6*abs(reference)`，physical limits 必須完全相同。最後一輪
須同時通過位移／速度門檻，先前輪次須未收斂並有合法 relaxation。

正常合成輸出與 15 種損壞注入已通過，包含缺少紀錄／rank、failed rank、
digest 損壞、rank history 分歧、錯誤 mapping、NaN、假收斂、超大 iteration
count 與數值／policy 差異，見 `strong-fsi-checker-negative-audit.json`。
這是檢查器協定驗證；實際 1／2／4-rank numerical matrix 仍待完成。

```bash
python3 scripts/hpc_check_strong_fsi.py \
  --reference-run outputs/hpc03/moving-graph-v1/strong-fsi-history-1 \
  --run outputs/hpc03/moving-graph-v1/strong-fsi-history-1 \
  --run outputs/hpc03/moving-graph-v1/strong-fsi-history-2 \
  --run outputs/hpc03/moving-graph-v1/strong-fsi-history-4
```

膜數值 payload 元件 `MembraneCheckpoint.hpp` 已加入。保存 model／state
identity、authority 提供的 accepted step／time、穩定 node IDs 與 scalar
位移／速度；讀回須符合外部預期的 model、state identity 與 clock，並經
模型自身的 shape／finite／clamp 驗證。只回傳 restore candidate，不直接
替換已接受狀態。預設 4096 nodes，受既有 metadata byte cap 約束；時鐘由
outer runtime／bundle 提供，core 本身不具 accepted-clock authority。

`make -C solvers/cpu membrane_checkpoint_test pretensioned_membrane_test` 已
無警告建置，兩個 executable 皆 exit 0。非零膜狀態逐 byte 往返一致，
更換 dt 與載荷後的下一步位移／速度／committed state identity 與不中斷
運行精確相同；錯誤 model／state／clock、截斷、trailing bytes、超限與
payload 損壞均拒絕且不改動模型。證據綁定於
`membrane-checkpoint-source-manifest.json`。這只驗證 numerical core 的
續算；single-owner wrapper clocks／publications、分散 fluid shards 與
atomic paired restore 尚未接線，不能視為完整 FSI checkpoint。

追加 stable node ID 損壞、clamped node 非零位移與正無限大 payload
拒絕案例，重新建置無警告、測試 exit 0；來源、binary 與 logs 綁定於
`membrane-checkpoint-guards-source-manifest.json`，前版來源與證據保留。

`SingleOwnerMembraneRuntime` 現提供 collective `CaptureCheckpoint` 與
fresh-runtime-only `RestoreCheckpoint`。僅 numerical owner 持有 payload；
外部 bundle 的 metadata SHA-256 綁定 configuration、numerical bytes、
accepted context 與原 publication producer。還原先驗證並建構 numerical
candidate，再依目標 partition 分送 kinematics，全部成功後才以 noexcept
swap 接受。獨立保存 committed context，避免 abort 後誤用 trial clock。

`fluid_surface_traction_test` 的 pressure／viscous 非零案例已在 1／2／4
ranks 通過，7 份 rank reports 均 exit 0、無 timeout，來源及 log hashes
已核對於 `wrapper-checkpoint-matrix-audit.json`。涵蓋 metadata／numerical
損壞及截斷的 collective 拒絕、健康 retry、publication identity 與 payload
精確往返、下一步 trial identity 與不中斷運行精確一致、active-trial capture
與重複 restore 拒絕，以及 abort 後保存最後 accepted clock。每次讀寫使用
相同 rank 數；跨 rank 數重分配、fluid shards 與 paired bundle 原子還原
仍待驗證／整合，本結果不代表完整 FSI restart 已完成。

後續跨作業驗證已完成 1→2、1→4、4→1 ranks：獨立 MPI processes
寫入／讀取 pressure 與 viscous 非零案例的膜 checkpoint，來源 numerical
owner 為最後一個 rank，讀取端改為 rank 0。5 次作業、12 份 rank reports
均通過；metadata、numerical payload 重新保存後逐 byte 相同，accepted
位移／速度及下一步結果與目標 rank 數的不中斷 runtime 比較滿足
`abs(a-b) <= 1e-12 * max(1,abs(a),abs(b))`。
`wrapper-repartition-matrix-audit.json` 綁定來源、checkpoint 檔案及全部
rank logs，`wrapper-repartition-matrix-runs.json` 保留完整執行命令。
測試使用 `-membrane_checkpoint_write DIR`／`-membrane_checkpoint_read DIR`
且 DIR 由 driver 預先建立；這些 fixture 檔案讀寫不是 production bundle
的原子發布協定。膜 wrapper 的跨 rank 數續算已驗證，fluid shards 與
paired bundle restore 仍待完成。

再於 restore candidate 已完整分送、尚未 swap 接受前注入 rank 0 失敗。
1／2／4 ranks 均 collective 拒絕，所有 rank 仍無 committed publication、
也不能 capture checkpoint；移除故障後可健康重試，publication 與下一步
trial identity 精確一致。相同執行亦讀取前版 1-rank checkpoint 驗證相容
續算。7 份 reports、來源與 logs 已核對於
`wrapper-restore-failure-matrix-audit.json`。故障入口僅在測試 macro 下提供。

非零 strong FSI 的 1-rank exhaustion／retry 長測試已完成：先限制為
1 次迭代並確認明確的未收斂錯誤、accepted state exact rollback，再用
正常控制完成兩步（6／7 iterations）。step 1 的 RMS／velocity×dt 為
`1.9011931674423223e-8`／`3.8023863348846446e-8`，threshold
`6.5417913784292356e-8`；step 2 為 `3.6118189854073674e-9`／
`4.5247501270868059e-8`，threshold `6.6069268614094811e-8`。
rank report exit 0、無 timeout，wall 2539.1984 s；此為功能測試時間，
不是效能基準。`strong-fsi-exhaustion-1-audit.json` 核對了 logs、兩步門檻
與 frozen binary 的所有來源（含修改前封存版本）。

新增 `MovingConservationCheckpoint.hpp` 保存完整已接受 transition 的
source／target geometry 與 publication IDs、clock、體積、全部 endpoint
aggregates／boundary-label maps、roundoff 計數及正規化診斷。解析須符合
外部 payload SHA-256 與 target geometry／publication／clock，且受既有
metadata cap 限制。實際 moving runtime 回歸已加入 byte roundtrip、錯誤
authority／geometry／step、截斷及 trailing bytes 拒絕，建置完成；MPI
matrix 正在執行，尚未據此宣告完整 fluid checkpoint 或 FSI restart 完成。

首次 conservation matrix 在 1-rank serial oracle KSP 失敗，尚未到達
checkpoint 檢查；`conservation-checkpoint-oracle-failure-audit.json` 保留
非零退出與 log hashes。後續使用 `deferred-abort-4` 已驗證的明確 oracle
MUMPS／GMRES、MUMPS ICNTL(14)=100 設定；數值驗收門檻不變。

`ImmersedMovingTransientDistributedRuntime` 已加入記憶體中的
`CaptureAcceptedCheckpoint`／`RestoreAcceptedCheckpoint`：owned field、
geometry／publication／active-layout／material identities、accepted clock
與完整 conservation record 一起保存。fresh target 必須先重建同一幾何
publication 與 index；驗證目標與 owned shape、配置候選資料後，經數值
setter 完成 PETSc staging，再接受守恆紀錄。已接受過 transition 的 target
不能重新 restore。這是 outer reader 使用的 typed candidate 介面，尚不
負責檔案認證、configuration binding、跨 rank shard 重分配或 paired commit。

回歸新增重建前後幾何、錯誤 layout／owned range／clock 的 collective
拒絕及 retry、owned field 精確往返、conservation byte roundtrip，以及
下一步求解與不中斷運行的 scaled L2 `<1e-8` 比較。`fluid-checkpoint-test`
及 `fluid-checkpoint-source-manifest.json` 已 frozen，1／2／4-rank matrix
正在執行，結果尚待核對。

restore preflight 再補強：完整 conservation metadata 先經 codec 驗證，
再 collective 比對 SHA-256，避免各 rank 接受不同的診斷紀錄；非有限
owned field／診斷、錯誤 source epoch 及 rank-local 紀錄差異均新增拒絕
案例。conservation type 移至獨立 header，讓 runtime 共用同一驗證邏輯，
並在任何數值向量發布前完成。新 build 無 compiler warnings，frozen
`fluid-restore-guards-test` 的 4-rank 故障／retry 測試正在執行；前版 matrix
仍使用其原 frozen executable，兩者證據不混用。

1-rank `fluid-checkpoint-test` 已通過：非零 uniform moving flow 的 owned
field 精確還原，守恆 metadata byte roundtrip，錯誤 layout／range／clock
拒絕後可 retry；下一步 scaled L2 為 0。`fluid-checkpoint-1-audit.json`
核對 rank report、logs 及所有 frozen 來源。這是 active layout 未改變的
記憶體 candidate 案例；2／4 ranks、額外 guards 與完整 bundle 尚待完成。

新增 `MovingAcceptedCheckpointMetadata.hpp`：shared metadata 綁定外部
configuration、geometry／publication／layout／material IDs、clock、global
row count 與 conservation SHA-256；解析只驗證目標並回傳 conservation
authority，不讀取 shards 或發布狀態。`moving_checkpoint_metadata_test`
已無 compiler warnings 建置並通過，涵蓋各 binding 拒絕、不同 owned
partition 的 shared bytes 相同、conservation authority，以及逐 byte 餵入
owned stream、截斷／trailing／錯誤 range／非有限值拒絕。證據綁定於
`moving-checkpoint-metadata-source-manifest.json`。

原 `BodyFittedAcceptedCheckpoint.hpp` 的三個 owned field stream 函式
原樣搬至 `OwnedCheckpointFieldStream.hpp`，介面與 wire bytes 不變。
既有 body-fitted accepted flow／transport checkpoint 回歸已通過三 ranks
及 split 1／2-rank communicators，含 steady／transient／outlet 與續算；
`shared-checkpoint-stream-audit.json` 核對全部 reports 與來源。

4-rank `fluid-restore-guards-test` 已通過：非有限 owned field／守恆值、
錯誤 source epoch 與 rank-local 守恆紀錄差異均在發布前 collective 拒絕；
healthy retry、owned field 精確還原及下一步續算成功，scaled L2 為
`2.45511e-15`。`fluid-restore-guards-4-audit.json` 核對四份 reports、logs
及 frozen 來源。本案例仍為 active layout 未改變的 typed-memory restore。

`OwnedCheckpointRepartitionReader` 可依任意順序串流讀取來源 shards，
只保留目標 owned rows，並驗證所有來源 global row ranges 的完整覆蓋、
無重疊及有限值。空 shard／zero-owned target 受支援；來源失敗後整個
reader 不可再接受結果。記憶體額外成本為 source range map，沒有 global
field buffer；目前每個目標 reader 仍掃描全部來源 bytes，I/O scalability
尚待評估，不能視為只讀取本地資料的效能路徑。

協定測試把 3 個非空 shards 加 1 個空 shard，重分配至 1／2／4／9
partitions，逆序逐 byte 輸入仍精確還原；亦拒絕重複、部分重疊、目標
owned range 外的缺口與截斷，失敗後 Finish 仍拒絕。build 無 warnings、
測試 exit 0，`owned-repartition-reader-source-manifest.json` 綁定來源及
binary/logs。這是 reader 協定驗證，尚不是實際 MPI 檔案 bundle 續跑。

原 `fluid-checkpoint-test` 的 1／2／4-rank matrix 全部完成，7 份 rank
reports exit 0、無 timeout；下一步 scaled L2 分別為 0、0、
`2.39334e-15`。`fluid-checkpoint-matrix-audit.json` 綁定每份 report、logs
及對應 frozen 來源。這仍是同 rank 數、active layout 不變的記憶體續算。

為重建保存時的 geometry publication，runtime 在 PrepareCommit 配置
上一個已接受材料表面的副本，FinalizeCommit 才接受；abort 清掉 trial
副本，保留 committed history。Capture 同時保存 previous／current
material，Restore 驗證材料身分、source clock 與目標幾何的 predecessor
material identity。多保留一份 replicated material metadata，不複製全域
fluid vector；既有 geometry／publication hash 算法未更改。

測試改由 material checkpoint bytes 往返後重建兩個幾何，新增 missing／
wrong predecessor 拒絕，再測下一步 active-layout 改變的續算。材料類別
不可 assignment，已改用 copy construction／emplace；最終 build 通過，
`checkpoint-material-history-test` 的 4-rank changing case 正在執行，結果
尚待核對。

`checkpoint-material-history-test` 的 4-rank changing case 已完成，四份
rank reports 均通過，owned field 精確還原；下一步 active layout 改變時
scaled L2 為 `1.85332e-14`，低於 `1e-8`。missing／wrong predecessor、
rank-local 失敗與 healthy retry 同時通過。完整來源及 logs 核對於
`checkpoint-material-history-4-audit.json`。

新增 `MovingCheckpointBundle.hpp`，沿用 immutable epoch／SHA-256／
receipts 協定，提供 domain shard catalog、local writer，以及先載入
previous／current material、再驗證 target metadata 並重分配 field shards
的 reader。shared files 包括 configuration／epoch／material binding、
accepted metadata、conservation 與兩份 material payload；流場各 rank
分片。這些 local I/O helper 不執行 MPI 或替 caller 發布 live runtime。

實際檔案協定測試已通過：來源 mapping `{2,0,1}`，三個 field shards
重分配至 1／2／4 partitions（含 zero-owned target）、材料與守恆往返、
錯誤 configuration／missing catalog／目標 owned range 外的檔案損壞拒絕，
修復測試檔案後可 retry。初次測試漏了 coordinator 的 receipt sorting，
被 Publisher 正確拒絕；補齊後以新目錄通過。最終 build 無 warnings，
來源、binary/logs 綁定於 `moving-bundle-source-manifest.json`。這是單行程
模擬來源 owners 的實際檔案協定測試，尚非 MPI numerical restart 或完整
FSI paired bundle；需要再接上 runtime factory 與 collective coordinator。

strong history matrix 的 1-rank 作業已完成，實際輸出通過
`hpc_check_strong_fsi.py` 的完整步／field／surface／history／port／守恆
協定檢查，保存於 `strong-fsi-history-1-self-check.json`。自我比較不提供
跨 rank 數值一致性的證據；2／4-rank 作業仍待完成。

1／2-rank strong history 比較已完成並通過：兩個 nonzero strong steps
的 17 組 fluid／surface field comparisons 最大 scaled L2 為
`1.2299955537051401e-12`，完整迭代歷史、ports 與 conservation 比較
最大 error/tolerance 為 `9.025606610871886e-4`。comparison 與 checker
hashes、所有 frozen 來源核對於 `strong-fsi-history-1-2-audit.json`；4-rank
作業仍待完成，不能提前視為全矩陣驗收。

新增 collective `RestoreMovingCheckpointRuntime`：先協議 manifest／
source mapping／configuration／path，載入材料並重建前後幾何，再配置
fresh runtime、串流讀取 target-owned field、驗證 predecessor geometry
與 conservation，最後接受候選。任何 field 階段錯誤都在各 rank 關閉
候選 runtime 後重拋；caller 仍須等其他 coupled candidates 成功才發布。

實際 MPI 回歸已接上每 rank shard writer、`GatherCheckpointReceipts`、
完整 catalog 發布、manifest discovery 與 factory 還原。新增 checksum
失敗（在候選 runtime 已配置後）、rank-local configuration 不一致，以及
原 accepted runtime 不變與健康重試。最終 build 通過，frozen
`file-runtime-factory-test` 的 4-rank changing case 正在執行；檔案 factory
與下一步 numerical continuation 的結果尚待核對。

檔案回歸增加 `-moving_checkpoint_read_only true` 與
`-moving_checkpoint_source_ranks N`，讀取端依保存時的 catalog／mapping
發現 manifest，並在新 communicator 上重分配場資料。來源 bundle 保持
唯讀，與目標 rank 數自行求得的 baseline 比較讀回場及下一步 scaled L2
`<1e-8`；同次 writer／reader 仍要求 owned field 與 conservation 精確往返。
checksum／configuration 拒絕與 healthy retry 保留。無 warnings 建置後
已 frozen `cross-rank-file-test`，4→2、4→1 的獨立 MPI 作業順序執行中；
來源 4-rank bundle 已完整發布，其原作業仍完成後續 continuation 驗證。

強耦合 1／2／4-rank history matrix 全部完成：7 份 rank reports 均 exit 0、
無 timeout，完整 field／surface／iteration／port／conservation 比較通過。
4-rank 對 1-rank 最大 field scaled L2 `1.5223305432149862e-12`，history
最大 error/tolerance `0.0013218419692356165`；全部來源、comparison、
checker 與 reports 綁定於 `strong-fsi-history-matrix-audit.json`。

`file-runtime-factory-test` 四份 reports 已通過，實際 shard 發布／discovery、
checksum 失敗後候選清理、configuration 拒絕與 retry、owned field 及
conservation 精確還原、changing-layout 下一步 scaled L2
`2.94854e-14` 均已核對於 `file-runtime-factory-4-audit.json`，包含實際
bundle 檔案 hashes。這份 frozen fixture 沒有 ports，不能據此宣告可變
port control 的 numerical restart 已驗證。

補上目前 queued port control values 的 capture／restore binding：metadata
payload `/2` 保存 bounded map；factory 在建立 fresh runtime 前驗證完整
port catalog 並套用保存值，immutable labels／control modes 由 case
configuration 保留。typed restore 要求新 runtime controls 與保存值相同，
並把 controls 納入 collective metadata agreement。`/1` 僅允許 port-free
cases 續跑；不能用初始化 control 值猜測舊 checkpoint 遺漏的資料。

新版 metadata／真實檔案 bundle 測試均通過，包含 controls readback、
套用保存值、錯誤／missing／nonfinite controls 拒絕且不改 candidate
options，以及 `/1` port-free 相容性與有 ports 時拒絕；證據綁定於
`port-checkpoint-source-manifest.json`。帶可變 ports 的完整 numerical
restart、跨 rank 數 file continuation 與 FSI paired restart 仍待驗收。

集中提交前以本階段完整來源重新建置 strong runtime，無 compiler
warnings；4-rank `strong-zero` 的 communicator／precommit／structure
prepare／invalid controls 拒絕、exact rollback 與兩步 healthy retry 均
通過，見 `checkpoint-stage-strong-audit.json`。這項 smoke regression
不取代前述 nonzero matrix，各 frozen 版本的證據仍分別保留。


跨 rank 數檔案矩陣已完成並核對：4→2-rank 的 owned restore／next-step
scaled L2 分別為 `1.33937e-14`／`3.44263e-14`；4→1-rank 分別為
`2.44193e-14`／`3.4676e-14`，均低於 `1e-8`。三份 rank reports
均正常退出且無 timeout；207 項來源身分、各 rank log hashes 與來源
bundle 檔案 hashes 綁定於 `cross-rank-file-matrix-audit.json`。
`owned_roundtrip_exact=1` 標記屬同次測試前段的 same-rank typed roundtrip，
跨 rank 的驗收依據是上述 scaled L2，不宣稱逐位元相同。這份 frozen
fixture 使用 port-free metadata `/1`，也不代表完整 FSI paired restart。

新增可選 `-moving_port_checkpoint_root` 數值回歸：接受非零 channel
FSI 步後把 inlet pressure control 從 `.1` 改為 `.2`，寫入真實 MPI
bundle，以原始 `.1` options 建立還原候選，確認保存的 `.2` controls
與既有 accepted field 均正確，並比較下一步流場。該 fixture 的 geometry
provider 合法地未附 predecessor；factory 因此先比較保存的 geometry／
publication identities，必要時重建不附 predecessor 的目標，再要求兩項
identity 完全一致。未更改既有 hash 定義。2-rank frozen executable 已啟動，
結果待核對；此測試尚未還原 structural runtime 或發布完整 FSI pair。


FSI adapter 新增 accepted publication checkpoint 介面，補足流場 owned
vector 與膜核心以外的 traction／consistent nodal force、producer／projection
stamp、composition identity 與 accepted coupling context。每份 payload
限定一個 owned surface slice，保存 stable node IDs；parser 以外部 payload
SHA、adapter configuration、材料／幾何身分與流場 accepted clock 驗證，
不把舊 partition 的 publication 直接換上新 stamp。格式限制全域 patch
節點數與 payload bytes，空 owned slice 仍有完整 stamp。

`RestoreCheckpoint` 只接受 fresh idle adapter，先驗證 fluid conservation
history，再解析各 rank 候選、協議共同 composition／材料／幾何／context，
最後才以 no-throw swap 發布。呼叫端仍必須等完整 fluid／structure pair
候選成功後一起發布；目前尚未加入 paired bundle 的 shard catalog 與
跨 rank surface publication 重建，不能視為完整 FSI 檔案續跑。

`moving_fsi_publication_checkpoint_test` 已通過精確 bytes／traction identity
往返，錯誤 configuration／geometry／material／interface／clock、checksum、
truncated／trailing bytes、nonfinite 值、節點上限、錯誤 partition 及合法
空 slice 測試。`fsi-publication-codec-result.json` 與
`fsi-publication-checkpoint-source-manifest.json` 保存來源／結果。無 compiler
warnings 建置後，2-rank `zero` adapter 回歸已啟動：測試 rank-local checksum
拒絕、解析成功後的跨 rank iteration 不一致拒絕、未發布狀態、healthy retry
與 overwrite 拒絕；數值程序結果仍待核對。


可變 port 的 2-rank 真實檔案數值回歸已完成：兩份 reports 正常退出、
無 timeout，約 519.9 秒。還原端採原始 `.1` inlet options，factory 正確
套用保存的 `.2` control；accepted owned field 精確保留，未附 predecessor
的 target publication identity 精確重建。原 runtime 與還原 runtime 的
下一步 scaled L2 為 `0`，還原端 commit 成功，原端 abort 後 accepted
field 不變。207 項來源身分、reports／log hashes 與實際 bundle 檔案
hashes 已核對於 `queued-port-checkpoint-2-audit.json`。這是同 rank 數
fluid candidate 驗收，未包含 membrane／adapter 的完整 paired restore。


2-rank adapter 回歸已完成並核對於
`fsi-publication-checkpoint-2-audit.json`：兩份 reports 均正常退出、無 timeout。
publication bytes／traction identity 精確還原，rank-local checksum 失敗與
已通過各 rank parser 後的 coupling iteration 不一致均拒絕，失敗後沒有
committed publication；healthy retry、overwrite 拒絕及原 fluid field
保持不變亦通過。既有 paired rollback／retry 與三個 decimal steps 均
通過。此測試使用 zero fixture 及相同 fluid owner，不代替 nonzero
paired file continuation；下一步仍是完整 paired catalog 與 fresh pair
候選建構、失敗清理及下一步強耦合驗證。


新增 `MovingFsiCheckpointBundle.hpp`，把 moving-flow catalog、每個 rank 的
FSI traction publication，以及單一 owner 的 membrane metadata／numerical
payload 放入同一 epoch。collective writer 核對 communicator、adapter 的
實際 fluid owner、來源 mapping 與完整 accepted coupling context；呼叫端
收齊 receipts 後發布完整 paired catalog。fresh-pair factory 先還原流場，
再建立 adapter／membrane、載入各 shard 並核對兩者 step／start／dt／
coupling iteration，全部成功才回傳 pair。後續失敗會 collective Close
候選流場，原 runtime 不受影響。surface publication 目前仍要求保存時的
rank count 與 ownership；流場的重分區能力不等同整個 FSI pair 重分區。

新增可選 `-moving_fsi_checkpoint_root` 強耦合回歸，在第一個 accepted
step 寫入完整 bundle，注入 membrane numerical checksum 錯誤以驗證
candidate cleanup／retry，再比較恢復的 field、adapter bytes 與 membrane
metadata／numerical bytes。第二步比較 fluid scaled L2 `<1e-8`、膜位移／
速度 absolute error `<1e-12`，traction／force、完整 Aitken history、ports
與五項守恆量使用 `1e-12 + 1e-6*abs(reference)`；迭代次數也必須相同。

第一輪 2-rank zero 回歸在第一個 strong step 成功後，因 writer 合併的
local receipts 未排序，被 `GatherCheckpointReceipts` 正確拒絕；manifest
尚未發布，不能算 numerical restart 通過。失敗來源／reports／log hashes
保存在 `paired-checkpoint-receipt-order-failure-audit.json`。已修正合併後的
排序，正在重建並準備以新的輸出目錄重跑，未改任何數值門檻。


修正排序後，完整 paired file 的 2-rank zero 回歸已通過。流場、adapter
publication bytes 與膜 metadata／numerical bytes 精確還原；膜 numerical
checksum 失敗後的候選清理及健康重試通過。下一個 strong step 一次迭代
收斂，原程序與 fresh pair 的 field scaled L2 為 `0`，surface／history／
ports／守恆比較全部通過。所有 rank 狀態、來源與 bundle／log hashes
已綁定於 `paired-checkpoint-sorted-zero-2-audit.json`。同一 frozen binary
的 2-rank nonzero strong case 已啟動，輸出在
`paired-checkpoint-nonzero-2`，結果待核對；兩項測試都在同一 MPI job 內
建立 fresh pair，不等同獨立作業或跨 rank 數的 paired restart。


成對 file fixture 增加 `-moving_fsi_checkpoint_read_only true`：新 MPI 作業
只 discovery／read 已發布的完整 bundle，不重新寫入來源；仍自行計算
基準軌跡，用於 accepted field 的 scaled L2 `<1e-8` 與下一步完整強耦合
比較。adapter publication 與 membrane metadata／numerical 的 recapture
直接對照保存的 shard bytes，要求精確相等。同次 writer／reader 另外
保留與原 runtime 全部 bytes／field 精確相同的要求，不用跨程序 bitwise
一致假設取代數值比較。

無 compiler warnings 建置後，`paired-checkpoint-readonly-test` 已凍結。
2-rank zero 新作業讀取前次已終止的 writer bundle，來源檔案的執行前
hashes 保存於 `paired-checkpoint-readonly-zero-input-hashes.json`，結果
待核對；非零 writer／fresh-pair continuation 作業仍在執行。


獨立 2-rank zero reader 作業已完成：accepted field 與下一步 field scaled
L2 均為 `0`，adapter／membrane 保存 bytes 精確 recapture，checksum
失敗後清理／重試與下一步 surface／完整 history／ports／守恆比較通過。
來源 bundle 的 12 個檔案集合與 hashes 前後完全一致；來源、writer
證據與 reader reports／logs 綁定於
`paired-checkpoint-readonly-zero-2-audit.json`。這證明 writer 終止後的新
MPI 作業可還原完整 zero pair，非零獨立作業仍待驗證。

非零 paired checkpoint 驗證擴展至 1／2／4 ranks：2-rank 使用
`paired-checkpoint-sorted-test`，1／4-rank 使用新增 read-only fixture 的
`paired-checkpoint-readonly-test`（本次兩者皆以 writer 模式執行）；物理
runtime headers 相同，各自 source manifest 保留。輸出目錄為
`paired-checkpoint-nonzero-{1,2,4}`，三組都仍執行中。分別配置本機
CPU 集合 2,3／0,1／4,5,6,7，OMP／OpenBLAS threads 均為 1；此並行
安排用於功能驗收，不用其 wall time 宣稱 strong scaling。


非零 paired file 的 1／2／4-rank 矩陣已全部完成，7 份 rank reports
正常退出且無 timeout。下一步均為 7 次強耦合迭代，field scaled L2
分別為 `0`／`0`／`7.6475960009046241e-13`，全部低於 `1e-8`；
surface、完整 history、ports 與五項守恆比較通過。完整 catalog 發布、
membrane checksum 失敗後 flow candidate 清理／retry，以及 accepted
flow／adapter／membrane 精確還原均通過。

原基準軌跡另外用 `hpc_check_strong_fsi.py` 比較跨 rank fields／history／
ports／守恆：2／4-rank 對 1-rank 的最大 field scaled L2 分別為
`1.2299955537051401e-12`／`1.7302690543569349e-12`。來源／archives、
checker、comparison、全部 reports／logs 與 bundle hashes 綁定於
`paired-checkpoint-nonzero-matrix-audit.json`。這仍是同 job 的 fresh pair
驗收；writer 已終止後的新 2-rank nonzero reader 已啟動，輸出在
`paired-checkpoint-readonly-nonzero-2`；結果如下。

上述非零獨立 reader 已完成：兩份 rank reports 正常退出且無 timeout，
accepted field 與下一步 field scaled L2 均為 `0`，保存的 adapter／membrane
bytes 精確 recapture，checksum cleanup／retry 及下一步 7-iteration
surface／history／ports／守恆比較通過。來源 bundle 的 12 個檔案 hashes
前後完全一致；來源、writer matrix 與 reader reports／logs 綁定於
`paired-checkpoint-readonly-nonzero-2-audit.json`。因此同 rank 數完整
nonzero pair 已證明能由 writer 終止後的新 MPI 作業續跑。

完整 pair 的 surface restore 現在可重分區：讀取每個已驗證的來源
publication shard，檢查共同 reference／來源 layout／composition／context，
並要求 stable node IDs 全域剛好覆蓋一次。之後只保留目標 rank 擁有的
traction／nodal force，以目標 partition identity 和全部來源 payload
identities 建立新的 producer／projection provenance；舊 partition stamp
不會被套到新 ownership。同 rank 且 owned IDs／layout／partition identity
完全相同時仍使用原 payload，保留精確 bytes 還原。fresh membrane owner
可與來源 owner 不同，既有 owner checkpoint restore 會把完整核心狀態
重新發布至目標 surface partition。

測試 CLI 新增 `-moving_fsi_checkpoint_source_ranks N`，read-only 作業以
來源 rank 數建立 compatibility／catalog，目標 MPI communicator 則自行
建立 surface ownership；accepted surface 及下一步的 field、surface、
history、ports、守恆均與目標 rank 基準比較。已用目前 frozen 4-rank
nonzero bundle 啟動 4→2-rank 驗證；來源 bundle hashes 已在執行前保存。

目前 repartition 來源另完成 2-rank zero smoke：publication exact／checksum／
typed-state／overwrite guards、paired prepare／context rollback、三個 decimal
steps，以及 strong communicator／precommit／structure／invalid-control 失敗後
兩步 retry 均通過。211 項 runtime 來源身分與四份 rank reports／logs 已核對
於 `paired-repartition-head-smoke-audit.json`；nonzero 重分區仍依 frozen
4→2／4→1 與後續 1→4 作業的獨立證據驗收。

首輪 frozen 4→2 作業已完成 checkpoint restore，但測試在續算前誤要求
來源 4-rank 與目標 2-rank baseline 的 active-layout SHA 逐位元相同而退出。
該 SHA 包含 geometry identity；不同 rank 的獨立非零求解可有已通過既定
數值門檻的 roundoff，因此不能拿此 hash 取代 stable active row catalog
檢查。兩份 failed reports、logs 與精確來源已保存在
`paired-repartition-layout-hash-failure-audit.json`，不列為重分區通過。
測試改為跨 rank 精確比較 active node IDs、port IDs 與 gauge presence；
同 rank restore 仍要求完整 layout SHA 相同。已停止尚未抵達同一錯誤檢查的
舊 4→1 作業，待修正版 binary 一併重跑；數值門檻未調整。

最新 publication codec guards 另以 ASan／UBSan 建置並通過，來源與
sanitized binary 綁定於 `paired-repartition-codec-sanitizer-audit.json`。
受管理環境的 ptrace 使 LeakSanitizer 在測試開始前無法啟動；因此本次
address／undefined-behavior 檢查明確使用 `detect_leaks=0`，不宣稱 leak
掃描結果。修正版 4→2、4→1 與反向 1→4 nonzero 作業目前並行執行。

新增 `scripts/hpc_check_paired_fsi_restart.py` 作為重分區驗收器：它會重新
驗證每份 rank report、stdout／stderr hashes、兩步完整 strong history、
checkpoint／continuation marker、來源／目標 rank 模式與輸入 bundle 前後
hashes，並拒絕超過 `1e-8` 的 accepted／next-field scaled L2。使用先前已
稽核的 nonzero logs 建立暫存相容 fixture，健康案例通過；把單一 accepted
error 改為 `1` 後確實非零退出，結果見
`paired-restart-checker-self-test.json`。該 synthetic self-test 只驗證 checker，
不新增數值證據；正式結論仍取自修正版實際作業。

修正版 nonzero matrix 最終全部通過。4→2、4→1、1→4 的 accepted field
scaled L2 分別為 `5.370410638499683e-13`、`3.303457280884551e-13`、
`5.958747452066621e-13`；下一步分別為 `1.431711374254415e-12`、
`1.679465531835937e-12`、`1.159705383816217e-12`，均低於原定 `1e-8`。
三組下一步都以 7 次強耦合迭代收斂，accepted surface、下一步 surface、
完整 history、ports 與五項守恆量全部通過。Checker 另重新驗證每份 rank
report、timeout、stdout／stderr digest 及來源 bundle 前後 hashes。

第一輪 1→4 暴露來源 layout 驗證錯誤：程式把目標 rank 的局部 ownership
切片拿來重建單 rank 來源 layout，因此四 rank 一致在 restore payload 階段
拒絕。該失敗保存於 `paired-repartition-fixed-1-to-4-failure-audit.json`，不作
數值通過宣稱。修正後由已認證來源 shards 彼此核對共同 layout identity，並
另用 reference identity 與 stable node ID 精確全域覆蓋綁定目標；不依賴不同
rank 分割下可能含無害 roundoff 的重建 geometry hash。修正版 1→4 使用
`paired-repartition-layout-validation-source-manifest.json`；4→2／4→1 使用
先前 frozen `paired-repartition-fixed-source-manifest.json`，最後修正不改變其
正常來源覆蓋與 provenance 路徑。

三份正式結果為 `paired-repartition-fixed-4-to-2-audit.json`、
`paired-repartition-fixed-4-to-1-audit.json` 與
`paired-repartition-layout-validation-1-to-4-audit.json`。完整命令、環境、RSS、
來源對應與限制整理於
[HPC-05D 完成報告](HPC_05D_MOVING_FSI_RESTART_REPORT.md)。這些作業與其他
工作負載並行，只提供功能驗收，不作 scaling 或跨節點結論。最終來源與
證據索引為 `paired-repartition-completion-audit.json`。
