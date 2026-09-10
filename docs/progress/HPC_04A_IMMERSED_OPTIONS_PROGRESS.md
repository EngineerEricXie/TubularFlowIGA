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
因此此 gate failure 在 options 改動前即存在，仍須定位修復，不能列為通過。
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
現有 `rigid` failure 尚未修復，HPC-04A 不勾選完成。

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

## 剩餘工作

已完成上述 81 個作業（64 正向、17 預期負向）與 131 份成功 rank report。
另有上述完整 FSI 回歸通過；`moving-regression-final-v2` 已終結且失敗，須定位並修復
rigid wall trace 誤差，再重驗及補上
整體驗收與 source/evidence 彙整；不得把正在執行的測試列為通過。
新增巢狀 viewer 與本機後端能力驗證見 [診斷進度](HPC_04A_BACKEND_DIAGNOSTICS_PROGRESS.md)。
舊版 `iga_transport` options 隔離已完成，見 [legacy options 驗證](HPC_04A_LEGACY_OPTIONS_PROGRESS.md)。HPC-04B／C 的其餘預條件器與大小案例評估仍待完成。
