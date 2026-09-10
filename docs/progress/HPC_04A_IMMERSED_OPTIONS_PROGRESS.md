# HPC-04A immersed solver options 開發紀錄

- 狀態：實作與驗收進行中，HPC-04A 尚未完成。
- 基準 revision：`e332973f5693ef4ae0a19548a950abbbac278443`。
- 範圍：serial／distributed static、fixed transient、moving epochs、FSI wrapper 與 native graph。
- 共用 family options 層已獨立提交：`bc5de4f5ae99b63afa565b5a41b1ddb942311e53`，
  見 [family options 紀錄](HPC_04A_FAMILY_OPTIONS_PROGRESS.md)。
- 操作介面：[SOLVER_OPTIONS.md](../SOLVER_OPTIONS.md)。

## 實作

Static／transient options 的尾端新增 `solver_options_prefix`，保留既有 aggregate
成員位置。空 prefix 沿用 `immersed_static_`／`immersed_transient_`；native graph 及
FSI 依 domain ID 產生 prefix。Domain 選項覆蓋 family 選項；immersed 延續不繼承
未加 prefix 的 root KSP／PC 設定。Serial transient 新增 family options 支援。

各 runtime 維持 private options snapshot，SetFromOptions／solve 的 options 與
error-handler scope 返回後恢復。Moving 的 committed／trial epochs 共用 immutable
owner。KSP 先銷毀，再釋放 options；幾何仍保留到對應 runtime 銷毀。
`SolverConfiguration()` 為本地查詢；native 每個 accepted step 輸出有效 KSP／PC。

共用 helper 增加可選 family fallback 與 root inheritance policy，並依 PETSc
不分大小寫的規則處理 key 優先順序；來源 value 與 usage flag 對應仍保留。
Transient input hash 包含正規化的 scoped options，版本更新為 fixed v6／moving v3。
沒有修改 mesh、database、場輸出、port 或 checkpoint payload 格式。

## 驗收方式

本機 GCC 11.4、OpenMPI 4.1.2、PETSc 3.15.5 real64／int32／MUMPS。
CPU tests 使用 C++17／O3／warnings；native graph 沿用 coupling Makefile 的 OpenMP
建置，執行時 `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`。
小案例在本機執行；沒有跨節點或大型效能宣告。

```bash
make -C solvers/cpu -j2 petsc_solver_options_test \
  immersed_distributed_static_flow_test immersed_transient_distributed_runtime_test \
  moving_immersed_solver_options_test moving_immersed_transient_flow_test \
  moving_immersed_transient_flow_fsi_runtime_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
make -C solvers/coupling -j2 iga_multidomain_flow immersed_case_factory_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
python3 scripts/hpc_immersed_solver_prefixes.py --suite core --output-dir NEW_CORE
python3 scripts/hpc_immersed_solver_prefixes.py --suite native --output-dir NEW_NATIVE
```

Harness 對每個 rank 保存 exit、timeout、RSS 與 stdout／stderr。Core fixture 既有的
serial／MPI 場、history、conservation、failure／rollback／retry 門檻保留。
Native accepted ports 使用 `1e-12 + 1e-6 * abs(reference)`；實際 configuration
必須符合 family/domain override，且最後求解 reason 為正。Moving solver 比較的
場 relative L2 門檻為 `1e-8`，不因失敗而放寬。

## 中間驗證與失敗紀錄

證據根目錄：`outputs/hpc04/immersed/`。以下屬中間版本，不能代替最後 source 的驗收。

- `core-v4` 的 static default／override 與 transient default／override／split 均通過。
  Static 場 relative L2 約 `3.34e-11`；transient 約 `1.22e-14`。
- `native-v2` 的 6 個 static／transient graph 作業通過，10 份 rank report。
- `negative-v1` 的 8 個無效 KSP／backend 作業均 exit 1、沒有發布結果。
- `reference-v1` 比較保存的修改前 native binary：static／transient、1／2 ranks 的
  32 份 accepted output 檔案逐位元相同。保存的 reference binary 位於
  `outputs/hpc04/one-d/native-v1`。
- `one-d-regression` 的四種 implicit 方法、override、checkpoint save／resume
  共 20 個作業、60 份 rank report 通過。
- 最初 core harness 把 PETSc flags 傳給 fixture 的 positional parser，被正常拒絕；
  改用 `PETSC_OPTIONS`。沒有更改 production CLI 介面。
- 最初 native 診斷解參考已移走的 `unique_ptr`，發生 SIGSEGV；已改用既有
  `immersed_audit` 查詢表，修正後的 native 作業通過。
- 新 moving 比較最初把 fixture 的 KSP rtol 從原有 `1e-16` 改為 `1e-11`，
  雖回報收斂，場誤差仍超過門檻。正在以原 tolerance 及明確誤差輸出重測。
- 單獨建置新 moving target 時，因 PETSC_TARGETS 宣告太晚而缺少 PETSc flags；
  已移至 Makefile 前段。`build-moving-v2.log` 的建置失敗，其後舊 binary 的結果
  不計為新 source 驗收。

## 最終來源已完成的驗證

Native source identity：`1ccfb5e9f336fd7e27b2dc699c8e7f08f53630b178f9b3f58ca87192d93707a6`。
保存的 binary：`outputs/hpc04/immersed/native-v4-binary`。以下通過項目使用最終 runtime
headers；moving 比較 fixture 的精度調整與其獨立結果另列，不能把先前失敗列作通過。

| 驗證 | 作業 | 結果／證據 |
|---|---:|---|
| 共用 options world／split | 1 | `unit-v4`：3 rank reports，42 次 solve，最大 L2 `5.44e-16` |
| Static／transient default、override、split | 5 | `distributed-final-v2` 的五個案例全部通過 |
| Native family／domain override | 6 | `native-v4`，1／2 ranks accepted ports 與 solver 設定通過 |
| 預設相容性及無效 solver／backend | 20 | `compatibility-final-v2`：12 正向、8 負向；32 accepted output 檔案逐位元相同，錯誤後健康 process retry 通過 |
| 四種 1D 方法／checkpoint | 20 | `one-d-final-v2`：60 份 rank report，default／override／save／resume 通過 |
| Standalone flow／transport | 27 | `standalone-final`：18 正向、9 負向，場與 checkpoint 比較通過 |
| Moving solver 比較 | 1 | `moving-final-v2`：兩步與跨 epoch snapshot 通過；最大場 relative L2 `2.68e-12` |
| ASan／UBSan aligned 模式 | 1 | `asan-final-v2`：兩步通過，最大 relative L2 `7.41e-15`，沒有 sanitizer 診斷 |

已完成作業的 source／binary 與 rank 統計先保存於 `interim-acceptance.json`，
狀態為 `partial_complete`；兩個完整回歸結束前不生成整體通過紀錄。
例如最終 static default rank 0 的診斷為 assembly `39.884 s`、linear solve
`0.306 s`；transient default 為 `1.498 s`／`0.144 s`。一般作業最大單 rank RSS
`110,948,352 bytes`，ASan 為 `638,824,448 bytes`。部分驗證作業同時執行，這些數字
只供功能回歸觀測，不作效能或 scaling 比較。CUDA、跨節點通訊與大型 I/O：N/A。

Moving 的原 KSP tolerance 重測揭示非線性停止條件是比較精度的限制：第一步 relative
L2 `5.03e-10` 通過；第二步為 `4.10e-8`，其中一個求解器的 nonlinear residual
`3.26e-12` 已符合原 absolute tolerance。Aligned sanitizer fixture 也出現
`1.03e-8` 場誤差與 `2.53e-12` residual。這些數值失敗不是 sanitizer 記憶體錯誤，
也不是通過紀錄。最終 solver 比較 fixture 改用 nonlinear absolute `1e-14`、relative
`1e-11`；保留原 `1e-8` 場門檻，production defaults 不變。另提供
`-test_solver_aligned`，讓同一測試在較小的 aligned geometry 下驗證 owner 生命週期。
最終兩步的切割與 aligned 比較均通過，確認停止條件收緊後能達到原場誤差門檻。

ASan 使用 `-O1 -fsanitize=address,undefined -fno-omit-frame-pointer`，執行時
`ASAN_OPTIONS=detect_leaks=0 UBSAN_OPTIONS=halt_on_error=1`。不宣稱已檢查第三方
MPI／PETSc 的 process-exit leaks。指令向量保存於 `asan-final-command.json`；執行
`moving-options-final-asan -test_solver_aligned`。早期大切割 sanitizer 作業 `asan-v1`
達 1200 秒 timeout，未計通過。

完整 FSI 的 `fsi-final-v1` 達 1200 秒 timeout，沒有數值錯誤輸出，未計通過；
`fsi-final-v2` 使用 3600 秒上限與 `-domain_fluid_flow_ksp_converged_reason` 重跑，
已於 1428.025 秒退出 0，`converged=true`，host peak RSS 73,764,864 bytes。
此為 revision `01aa4c8` 的完整既有 FSI fixture，原驗收條件保留。
完整 moving 的 `moving-regression-final-v1` 達 1800 秒 timeout，未計通過；
`moving-regression-final-v2` 保留原 fixture，改用 7200 秒上限與
`-immersed_transient_ksp_converged_reason`，於 4924.678 秒退出 1（非 timeout），
失敗點為 `rigid target quadrature u-w gate failed`。Peak RSS 380,329,984 bytes。
獨立抽出的同一 rigid fixture 已比較 `bc5de4f`（改動前）與 `01aa4c8`
（改動後 runtime）：兩者皆退出 1，`u-w` gap 完全相同，為
`8.998129639693216e-6`，原門檻為 `5.684341886080802e-14`。
兩者 nonlinear residual 均為 `1.0010909032568383e-14`，2 次 Newton、48 次 KSP。
因此此 gate failure 在 options 改動前即存在；下述 fitting 修復已通過 focused rigid 驗收，完整回歸仍待重驗。
證據位於 `outputs/hpc04/rigid-debug/{before-v1,after-v1}/rank-0/`，
`source.json` 保存測試來源、兩個 binary 與編譯紀錄的 SHA-256。
測試共用原 fixture helpers、相同參數與環境，僅抽出 rigid case 並輸出數值；
原 `256*epsilon` gate 未放寬。

### Rigid failure 的積分定位

`rigid-debug/trace-v1` 保留原參數與 gate，結果為預期的既有失敗，非 timeout。
轉移後節點速度與原常數完全一致（最大誤差 0），求解前 surface `u-w`
gap 為 `5.5293885617803503e-17`；求解後升至上述 `8.998129639693216e-6`。
第一次 Newton residual 為 `6.2062113278804282e-5`，第二次為
`8.2403680970955447e-8`。因此誤差不是初始場轉移造成的。

另以 `quadrature-v1` 在同一 target geometry 直接累加每個全域 spline basis 的
`∫Ω grad(N_a)·w dV − ∫∂Ω N_a (w·n) dS`，使用 runtime 原體積及表面規則，
沒有 KSP／Newton 求解。缺陷 L2 為 `6.2062113278804906e-5`，
最大節點缺陷為 `1.8973701844825705e-5`，全域帶符號總和為
`7.297457341648011e-19`。L2 與初始 nonlinear residual 相符，
但總和幾乎為零：只驗證總流量不足以發現這個 basis-wise 積分不相容。
目前體積規則採 octree 內部採樣，表面規則則積分實際裁切面；
此組合在此 fixture 不滿足常數速度所需的離散分部積分恆等式。

`trace-source.json` 保存兩個診斷 source／binary／build log SHA-256 與 rank reports。
`quadrature-v1` 的退出 0 只表示量測完成，不代表此缺陷通過物理驗收。
下一步須修復積分一致性並驗證 conservative mixed form、Jacobian、既有靜態／移動
案例；不能以放寬 rigid gate、修改常數初始場或略過測試取代修復。
完整 moving 測試已補上高精度 seed／solved gap 與求解殘差輸出，保留所有原 gate。
其 binary 會因診斷重建而改變；前述 `01aa4c8` 完整失敗的執行身分仍以原 run 記錄為準。

### 積分修復的體積矩基礎

新增 `PolyhedralVolumeMoments.hpp`，由已驗證的 `ClosedTriangulatedSurface`
計算物理體積上的正規化 tensor monomial moments，逐軸最高六次。
使用 x 方向反導函數、散度定理與 triangle Duffy 映射；一維 Gauss 階數依
多項式總次數決定，並以 long double 補償加總降低帶符號邊界項的消去誤差。
原理參考 [Quadrature-free immersed isogeometric analysis](https://arxiv.org/abs/2107.09024)。
此實作處理封閉多面體，無凸性假設；它只提供目標矩，不產生帶負權重的 runtime 規則。

`make -C solvers/cpu polyhedral_volume_moments_test` 與
`solvers/cpu/polyhedral_volume_moments_test` 通過；1,029 個解析矩比較涵蓋所有
`0 <= a,b,c <= 6` 的四面體、非凸 L 柱體，以及大平移／非等向縮放／反向輸入面序。
最大相對誤差 `9.75782e-19`（本機 long double）；另驗證超出次數、零尺度及非有限原點拒絕。
測試使用解析體積積分，未以另一個相同演算法的結果作 oracle。

這是積分修復的基礎，尚未接入 cut-cell catalog：每個裁切 cell 的目標矩已由下述截斷反導函數提供，
仍須由目標矩求得符合既有有限正權重契約的規則，並驗證體積／表面分部積分、cap／失敗處理、
compact／expanded 路徑與 static／moving／distributed 結果。
此階段尚未修復 `rigid` failure；後續 catalog／rigid 結果見下節，HPC-04A 仍不勾選完成。

### Cartesian cell 的目標矩

新增 `PolyhedralBoxVolumeMoments.hpp`，計算封閉多面體與任意軸向 box
交集的正規化 monomial moments。令 `xi, eta, zeta` 為 box 正規化座標，
採用只有 x 分量的反導函數：

`F_x = dx/(a+1) * clamp(xi,0,1)^(a+1) * eta^b * zeta^c`

並將 y／z 支撐限制於 `[0,1]`。其散度是 box 內所需多項式、box 外為零。
在原封閉三角面上積分 `F_x n_x` 即得交集體積矩；x 大於 box 上界的面仍須保留，
其飽和反導函數提供截面項。因此不能只查詢與 cell AABB 相交的三角面。
實作裁切 y／z 支撐，再於 x 上界分段，避免非凸截面或多個截面環的 cap 三角化。
剛好落在 x 上界的完整 facet 只分配至一個分支，避免重複計數。

新增 686 個解析矩比較，涵蓋非凸 L 柱體與 box 交集、斜面四面體的 x-slab，
仍逐軸最高六次。最大相對誤差 `8.13152e-18`。另驗證 32 個相鄰 boxes 的
體積及 x 一階矩加總、完全內部 box（原面不與 cell 相交）、外部空 box、
重合邊界，以及 triangle／quadrature caps 與無效 bounds 拒絕。
原先 1,029 個全域矩測試亦通過。命令同上一節；本次 build／test logs 為
`outputs/hpc04/rigid-debug/box-moments-{build,test}.log`。

此路徑目前採 long double 裁切，不宣稱 exact predicates 或任意病態幾何下的
誤差保證；runtime 尚未使用它。接入前仍須處理精確裁切／退化幾何、正權重擬合、
目標矩殘差 gate 與完整失敗傳遞，並重驗原 rigid gate。

### 正權重擬合的實際 rigid cell 實驗

`PolyhedralBoxMomentOptions` 新增 box 正規化座標內的 `coordinate_origin`／
`coordinate_scale`，預設仍為原 `[0,1]^3` monomials。積分直接使用新 frame 的
反導函數，避免由舊 moments 做高條件數的事後座標轉換。新增 343 個置中／縮放解析矩
（含奇次與零矩）與零尺度拒絕；全域 1,029、box 1,029 個解析比較皆通過。

已匯出原 rigid target 的 26 個 cut cells、原 quadrature nodes 與六次 tensor moments：

- 原座標下的最小加權修正有 25 個 cells 產生非正權重，不能使用。
- 在節點支撐上置中／縮放後，11 個 cells 仍確實不足 343 秩；20 個 cells 有非正權重。
  不能將問題歸因於 solver tolerance，也不能只修正既有權重。
- 新 frame 的邊界目標矩與 cube 解析積分相比，最大絕對差 `1.62631e-17`、
  最大相對 L2 差 `6.97482e-17`。目標矩本身並非上述擬合失敗的來源。
- cell 23 的 8×8×8 支撐內候選點 NNLS 仍失敗（最大矩殘差 `1.75453e-4`）；
  擴充至 16×16×16、frame 座標乘 1.25，再排除 cell／物理域之外的點後通過。
- 相同擴充策略在全部 26 cells 通過；每個保留 343 個正權重點。最大矩絕對殘差
  `9.43690e-16`、最大相對 L2 殘差 `2.38848e-15`、最小保留權重 `1.17924e-8`。

將這些權重用於原 `EvaluateBasis` 與原 surface rules，直接重測
`∫ grad(N_a)·w − ∫ N_a(w·n)`：L2 降至 `4.68312e-18`、最大節點缺陷
`1.26335e-18`，通過 `256*epsilon` 獨立 gate；原值為 `6.20621e-5`。
這證明此實驗修正了實際 basis 的積分不相容，尚不代表完整 moving regression 通過。

證據為 `outputs/hpc04/rigid-debug/{moment-data-v2,augmented-all,quadrature-corrected-v1}`；
`positive-fit-audit.json` 保存來源、binary、矩／節點／權重與逐 cell 結果的 hashes。
所有失敗擬合保留於 `fit-v1.json`、`fit-v2.json`、`augmented-v1.json`。
NNLS 原型使用 NumPy 1.21.5／SciPy 1.8.0，未新增 runtime Python 依賴。
候選點的域內判斷在此實驗使用已知 cube 邊界；正式版本須使用既有驗證過的
surface predicates，並加入 C++ 受限擬合、資源 caps、矩殘差拒絕與完整回歸。

### C++ 非負擬合器

新增無外部線性代數依賴的 `NonnegativeLeastSquares.hpp`。以 column／RHS scaling、
兩次正交化 QR 插入、Givens QR 移除與 active-set 非負線搜尋求解，使用 long double
工作陣列，最後以原矩陣重新檢查實際回傳的 double 權重殘差。
`within_tolerance` 必須由呼叫端檢查；無法符合精確矩的輸入不會被誤報為通過。
明確限制 rows、columns、內部 workspace bytes 與 iterations，非有限輸入／不可表示
係數拒絕。工作記憶體上界不包含呼叫端已持有的矩陣與 RHS。

`make -C solvers/cpu nonnegative_least_squares_test` 與執行檔通過：
解析可行／不可行解、零 RHS、重複／零 columns、極端尺度、200 組獨立 KKT
最優性檢查、active-set 移除及各類 cap／輸入拒絕。
ASan／UBSan 同一測試通過且保持 leak detection；sandbox 版本在測試結束後因
LeakSanitizer 不支援 ptrace 而退出 1，該紀錄保留，未算通過；非 sandbox 重跑退出 0。

同一 26 個 rigid cells 的候選點與目標矩改由 C++ 擬合，全部通過 `rtol=1e-14`。
每個回傳 343 個正權重，獨立 long-double monomial 重算的最大相對矩殘差
`7.62447e-17`、最小正權重 `5.80472e-8`。最多 577 次 active-set 操作，
單 cell 最多 117 次 QR 移除；保守 workspace 上界最大 11,306,272 bytes。
26 cells 作業 wall 9.430 秒、peak RSS 19,001,344 bytes；不是跨實作速度比較。

以 C++ 權重再次量測原 basis 的連續方程積分缺陷，L2 為 `5.19439e-18`、
最大節點缺陷 `1.45045e-18`，通過獨立 `256*epsilon` gate。
證據為 `rigid-debug/{cpp-fit-v1,quadrature-cpp-v1,cpp-fit-independent.json}`，
`cpp-nnls-audit.json` 保存 source／binary／權重／test logs hashes 與 run reports。

正式 runtime 尚未接入：仍須以 surface predicates 產生合法候選點，
整合逐 cell moment fitting 與 cap／失敗處理、更新規則 diagnostics／hash 身分，
再重驗完整 rigid、static、moving 與 distributed 案例。

### C++ 逐 cell 規則建構器

新增 `FittedCutCellVolumeRule.hpp`，串接原 seed 的支撐 frame、Gauss 候選點、
`SurfaceSpatialIndex::LocatePoint`、邊界體積矩與 C++ NNLS。每個新節點必須明確
判定為 `Inside`；`Boundary` 不保留，`Ambiguous` 拒絕。候選階數可逐次增加，
耗盡後拋出失敗，不退回原本不相容的 seed 規則。
seed 數量、幾何查詢、候選點、矩陣配置、NNLS workspace／iterations 都有上限；
矩陣依最終候選數一次配置。回傳前以 long double 重算實際 double 節點／權重的矩殘差。

`make -C solvers/cpu fitted_cut_cell_volume_rule_test` 與執行檔通過：
非凸 L 柱體的全部 343 個解析矩、每個回傳節點的域內判斷與正權重、
query／seed／candidate caps、無效階數／空 seed，以及候選耗盡的拒絕。
最大解析矩絕對誤差 `1.44877e-17`。`make mesh-test` 亦通過。

原 rigid target 的 26 個 cut cells 已全部使用此建構器，沒有 Python 權重或
cube 特例 membership 判斷。總計 47,612 次 surface point queries；所有 cells
在第一個候選階數通過，各保留 343 個正權重，最大相對矩殘差 `7.36155e-17`。
以原 `EvaluateBasis`／surface rules 重測的連續方程積分缺陷 L2 為
`5.38696e-18`、最大節點缺陷 `1.49795e-18`，通過 `256*epsilon` 獨立 gate。
作業 wall 12.530 秒、peak RSS 53,657,600 bytes，包含原 geometry 建構與 basis audit。

證據在 `rigid-debug/quadrature-builder-v1`；`cpp-rule-builder-audit.json`
保存來源、binary、build／test／mesh logs hashes 與 run report。
此建構器尚未修改 catalog；剩餘工作是 expanded／compact 表示、診斷與 hash 身分整合，
再驗證原 rigid solve 與完整 static／moving／distributed 回歸。
邊界 moment clipping 仍是 long double 路徑，不宣稱 exact-predicate 積分誤差保證。

### Compact 正權重節點表示

`CompactCutCellVolumeRule` 新增獨立的 `fitted_points` 表示。非格點正權重規則
不能再編碼為原 octree block／sample mask；新表示直接保存小型擬合規則的節點與權重，
並禁止與 block／sample records 混用。既有 octree 建構器目前仍只產生原表示。

共同 point-count、iterator 與 validation 已支援此表示：要求有限、正權重、
`[0,1]^3` 座標與嚴格遞增節點順序，拒絕重複／逆序節點及混合表示。
非凸 L 案例的 343 個解析矩同時檢查 expanded 與 compact；另逐點比較座標與權重的
浮點位元以及遍歷順序，全部相同。無效表示的拒絕測試亦通過。
`make mesh-test` 通過；完整既有 `cut_cell_volume_quadrature_test` 已退出 0 並通過
（`cfe6ca0` 表示變更，binary/source 身分見 `compact-fitted-source.json`）。build／test logs 與 binary/source hashes 保存於
`outputs/hpc04/rigid-debug/compact-fitted-*` 與 `compact-legacy-test.log`。

尚未由 catalog 發布新表示：其 record／retained-byte diagnostics、hash 版本、
expanded／compact 共用的 fitting seed 與失敗傳遞須在接入時一併完成。
本次不宣稱完整 runtime 回歸已通過。

### Compact 計數與 hash 接線

共同 record 計數與容量計算已包含 fitted points，並檢查乘法／加法溢位。
`ValidateCompactStoredRule` 使用此計數，不能漏記新表示的 records 或 vector capacity。
容量數字依舊是 record backing storage，並非整個 process RSS。

`MovingCutGeometry` 改用共同 compact hash appender：原 block／sample 表示的 bytes
不變；fitted 表示以前置保留標記（合法舊 depth 不可能使用的 `0xffffffff`）與版本 1
區分，並綁定所有節點座標及權重。舊預設 geometry hash 不因空 fitted vector 而改變。
測試直接比較舊 block 及 sample 的原 byte stream，並驗證改動 fitted 座標或權重會改 hash。
record／容量計數及 overflow 拒絕、非凸解析矩與逐點一致性測試通過。

`moving_cut_geometry_test` 與 `make mesh-test` 通過；證據為
`rigid-debug/compact-hash-geometry-v1`、`compact-accounting-final-{build,test}.log`，
以及 `compact-accounting-audit.json` 中的 source／binary／log hashes。
Catalog 尚未建立 fitted 規則；下一步仍須加入共同 seed、規則發布與 fitting 診斷，
並重驗真正的 rigid flow solve。

### Catalog 發布與 rigid flow 修復

`CutCellVolumeQuadratureCatalog` 的第四個參數可提供 `FittedCutCellVolumeRuleOptions`。
正的 Cut cells 使用共同 compact octree seed，以串流方式擷取支撐 extrema；expanded
與 compact 不再各自選 seed，也不建立完整 expanded seed 陣列。發布前檢查原 octree
體積上下界、新規則的實際體積、正權重、output／logical／record caps 與 record workspace。
更新 published point counts、retained capacity、fit queries／candidates／iterations；
失敗發生在 constructor 的 cell 發布之前，不退回不符合矩的舊規則。

Moving geometry 以 `options.volume_fitting.emplace()` 選用此模式。
預設仍使用原 octree 模式與 v5 geometry hash；啟用 fitting 時使用 v6，綁定有效 fitting
參數、work counters 與規則內容。此選項尚未延伸至 native graph JSON 配置。

非凸 L catalog 測試同時建構 expanded／compact，逐點位元完全相同，物理體積通過解析值。
Compact 測試刻意設定 `max_points=1`，仍成功使用 streamed seed 與新表示；expanded
同上限會在發布前拒絕。測試另涵蓋 seed cap 與無效 fitting 配置拒絕。
此非凸 catalog 案例需要較大的候選集合，使用 `fit.max_columns=32768`；初次 8192-cap
失敗保留於 `catalog-fit-test.log`，未放寬矩殘差容差。

原 rigid fixture 已選用 moment fitting，其幾何、初始速度、時間步、KSP／Newton 容差
及 trace／conservation gates 保留。抽為共用 `RunRigidTranslation`，完整 moving 測試
與新 `moving_immersed_rigid_translation_test` 使用同一段驗收程式。
後者 expanded／compact 均通過原全部 trace 與 conservation gates：`u-w=5.52939e-17`，
門檻 `5.68434e-14`，nonlinear residual `3.97034e-16`，Newton／KSP iterations 均為 0。
這是正確常數初始場的離散殘差已消除；沒有改動驗收門檻或強制覆寫求解結果。

`fitted_cut_cell_volume_rule_test`、`moving_cut_geometry_test`、上述 focused rigid 與
`make mesh-test` 均通過。來源、binary、build／test logs 與 rank report 由
`rigid-debug/catalog-fitting-audit.json` 彙整。完整 moving executable 已重建，
接續作業為 `outputs/hpc04/immersed/moving-regression-fitted-v1`，尚未有完成結果。
六次 polynomial moments 的驗收不代表非線性／stabilization 積分完全精確；
其餘非剛體、static／distributed 比較與更廣泛的配置驗收仍須繼續。

### Geometry identity 與失敗建構相容性

新增 `moving_fitted_geometry_identity_test`，使用固定單 cell cube fixture，並以
pre-catalog revision `53b3238` 實際執行取得的 v5 digest 作為預設模式相容性基準：
`f2b2333e955112adb2663fb95214d01981b9fa2f195bd382b8809d44ee8140b1`。
現行預設模式維持相同 digest；啟用 fitting 後 identity 不同，重複建構則一致。
只增加未耗盡的 fitting iteration cap，規則座標／權重逐位元不變，identity 仍改變，
確認有效資源政策有納入身分。

無效 query cap 的下一次建構會拒絕；先前 geometry 的 identity、規則大小、
座標與權重均維持原值。測試以建構失敗前的 points 副本直接比較，未只依賴 cached hash。
此驗證涵蓋配置預檢失敗，不代表所有中途 fitting 失敗或 runtime rollback 已窮盡。
測試 exit 0；來源、binary 與 logs 見 `rigid-debug/fitting-identity-audit.json`。
完整 moving 回歸仍由獨立 frozen `c4fee00` executable 執行，尚未列為通過。

### Fitting 空規則與資源上限拒絕

補強 fitting constructor 的零支撐分支：只有 certified／unresolved volume 均為零的
seed 可以作為空規則返回。對仍有未解析體積的零取樣 seed，立即拒絕建構。
原預設 octree 模式仍保留其既有「建構後由 consumer validation 拒絕」行為。

`fitted_cut_cell_volume_rule_test` 新增 expanded／compact 的三類驗證：有效但耗盡的
query budget 與不足的 record workspace 拒絕；僅面接觸、物理體積確為零的 Cut cell
保持空規則且不執行 fitting；微小內部非凸幾何在 depth 0 無取樣支撐且 unresolved volume
為正時，fitting constructor 拒絕。既有 343 解析矩與兩種儲存模式逐點一致性亦通過。
首次邊界 fixture 未包含完整 surface，正確被背景域檢查拒絕，記於
`catalog-boundaries-test.log`；修正為包含 surface 的網格後，
`catalog-boundaries-v2-test.log` exit 0，未改動背景域驗證。
證據 hashes 見 `rigid-debug/catalog-boundaries-audit.json`。

### Fitting 發布的累積工作量

新規則發布時，現在把 fitted points 納入 `attempted_record_attempts`、
`attempted_logical_output_points`，以及 expanded 模式的 `attempted_output_points`。
上限涵蓋 compact seed 建構（含 rescue）與新規則發布的合計工作量；不再只檢查
最後一次 seed 的 records 或單獨的 fitted point count。此變更只影響可選 fitting 模式。
對累積 logical-point cap，測試設為 seed 與新規則總數減一；即使二者個別都低於上限，
expanded／compact 仍須拒絕。

更新後 catalog 全測試 exit 0；同一測試用 `10812e0` 的舊 catalog header 編譯，
exit 1 並明確報告 `fitted publication omitted cumulative records`，證明能抓到原漏計。
`moving_fitted_geometry_identity_test` 亦通過，舊預設 v5 digest 保持相同。
Fitting identity 綁定 diagnostics，故其新 digest 會反映正確累積計數；正在執行的完整
moving 回歸仍是 frozen `c4fee00` 的獨立證據。Logs 與 hashes 見
`rigid-debug/catalog-work-audit.json`，本次沒有更改 quadrature points 或權重。

## 完整 moving 回歸完成

`moving-regression-fitted-v1` 的 frozen `c4fee00` binary 已完成：exit 0、無 timeout，
wall 5582.585 s、peak RSS 356077568 bytes，stderr 為空，最終
`moving_immersed_transient_flow_tests=passed`。這是原完整測試，不是 focused rigid
替代：包含 contraction／expansion layout、moving conservation、狀態 rollback／
prepare／abort、stationary／mixed port capture 與 rigid translation 等原斷言。

Rigid wall trace gap=5.52938856178035e-17，原 tolerance=5.68434188608080e-14；
nonlinear residual=3.97033685409117e-16，Newton／KSP iterations=0。剛體初始場已
滿足離散方程，不需要以迭代修正。原 volume／wall／Reynolds／divergence theorem
門檻保持，僅此 fixture 明確選用前述正權重 fitting；其他 fixture 保留原積分選擇。
完整測試的 contraction rows 2919→1375，expansion 1375→2919。

`outputs/hpc04/immersed/moving-regression-fitted-v1/audit.json` 核對 binary SHA256、
stdout／stderr hashes、rank report、162 個來源檔案與 `git show c4fee00:PATH`
逐項相等。這不是目前 HEAD 的全回歸；後續 geometry identity／empty seed／累積
工作量修補由前節 focused tests 證明。舊 `moving-regression-final-v2` 的 rigid
失敗與其 frozen `01aa4c8` 證據保留，不改寫成通過。

此工作站 OMP／OpenBLAS=1，沿用 LU/MUMPS root options 與
`-immersed_transient_ksp_converged_reason`；family options 規則仍由各自測試覆蓋。
期間另有 duct 等回歸，這個 wall／RSS 不作無干擾效能比較。仍不是分散式 FSI
或跨節點的驗收證據。

## 剩餘工作

上述 81 個作業（64 正向、17 預期負向）、131 份成功 rank report，以及完整
FSI 回歸、修復後完整 moving 回歸均已有通過證據。HPC-04A 的整體完成判定仍須
合併所有 family／CLI／backend 的來源與範圍稽核，不能把各時期的測試視為同一
HEAD 完整執行。新增巢狀 viewer 與本機後端能力驗證見
[診斷進度](HPC_04A_BACKEND_DIAGNOSTICS_PROGRESS.md)。舊版 `iga_transport`
options 隔離見 [legacy options 驗證](HPC_04A_LEGACY_OPTIONS_PROGRESS.md)。
1024-element duct 的 LU／block-Jacobi 亦已通過原門檻，見
[較大方管進度](HPC_04B_LARGE_DUCT_PROGRESS.md)；HPC-04B／C 其他條件與
無干擾 mesh／rank／硬體擴展仍待完成。
