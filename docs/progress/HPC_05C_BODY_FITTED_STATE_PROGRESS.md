# HPC-05C：貼體 3D accepted-state 恢復介面

後續進度：本報告的磁碟 codec／新作業恢復缺口已由
[3D 分片驗收](HPC_05C_BODY_FITTED_BUNDLE_PROGRESS.md) 接續；下文保留本批原始範圍。

日期：2026-09-09。基準 revision：`49ae5c5eead287601bc0b3d0ff45766bfd2b8975`
加本批修改。狀態：**部分完成，HPC-05C 保持未勾選**。

本批補上貼體 flow／transport runtime 的記憶體狀態 capture／restore。恢復後的
三步計算在三 rank WORLD 與獨立的 1／2-rank 子群都與不中斷參考精確一致。
這批尚未加入 3D 磁碟 codec、bundle provider、graph CLI 或新 MPI 作業的載入；
不能將同作業內建立 fresh runtime 的比較稱為新程序 restart。

## 介面與狀態範圍

[OwnedCheckpointVector.hpp](../../solvers/cpu/include/OwnedCheckpointVector.hpp)
提供 local owned-row capture、validation 與 scratch staging：保存 global row count、
ownership interval、該 rank 的 FP64 值。驗證實際 Vec ownership／local size、形狀與
有限值；沒有 gather 完整場。操作本身不呼叫 collective，runtime 必須把它們放在
`CollectiveLocalStage` 裡，在進入下一個 collective 前一致傳播局部錯誤。

[TransientFlowRuntime](../../solvers/cpu/include/TransientFlowRuntime.hpp) 新增
`FlowAcceptedCheckpointState` 與 accepted capture／fresh restore，保存：

- owned `(u,v,w,p)`、完整 resolved boundaries、pressure tractions 與 outlet 動態值；
- accepted step count／time、graph macro dt、累積 linear iteration count；
- caller 提供並驗證的完整 configuration identity SHA-256。

accepted clock 只在 finalize 更新，abort 不會以 stale `trial_time_` 覆寫它。
checkpoint-bound owner 的下一步必須接續正確 step 與 time；防止 step index 與
累積 linear iteration count 的有號溢位。新介面拒絕 initial、active 或 closed owner，
只允許 fresh Committed candidate restore。shape／clock／model／boundary topology
先完成 local validation，然後核對 replicated state 的 SHA-256；最後才 staging、
重建當前場的 snapshot／BE history 並發布。新作業的 matrix／KSP／halo 仍由 runtime
建立，下一次 Begin／Solve 重建必要資料，不保存前一 trial 的 scratch history。

**graph macro dt 與 numerical dt 分開。** 穩態 flow kernel 的 numerical dt=0，
graph 仍可有正的 macro dt。constructor 最後的可選 checkpoint identity 與 macro dt
參數建立此綁定；transient owner 要求 macro dt 與 numerical dt 相同，steady owner
必須明確提供正的 macro dt。省略 identity 的現有 caller 保留舊 step／time 行為，
不能使用新 checkpoint API。第一次測試把兩個 dt 混用而被 kernel 拒絕，已修正並
保留失敗證據，沒有改成把穩態求解當作暫態。

[TransientTransportRuntime](../../solvers/cpu/include/TransientTransportRuntime.hpp)
新增 `TransportAcceptedCheckpointState`，保存 owned scalar vector、真實 accepted
step count 與 configuration identity。restore 恢復 current／committed fields 及兩個
step counters，下一步的 warm-start 決策與步數延續正確。`BeginStep` 拒絕 `INT_MAX`
計數繼續加一。原 `ReadState(path)` 設 steps=1 的 legacy 行為保留，新 provider 必須
使用 typed restore；沒有將舊 metadata 的不足悄悄帶入新介面。

以上 identity 是 **provider 的責任**：須涵蓋 domain／system、完整配置、外部輸入、
`.ntiga` 分區、node／field 順序及執行映射。本批不以 vector 長度相同替代身分檢查。
測試使用明確分開的 fixture identity，尚未實作 native graph 的 compatibility builder。

fresh runtimes 必須在 live graph 之外建立。所有 local validation 與群組 agreement
發生在 staging 前；若後續 PETSc staging／copy／swap 發生錯誤，provider 應丟棄整個
候選 owner。這不是對 live graph 多個 PETSc Vec 的無條件原子替換承諾；全作業
candidate agreement 與 owner publication 仍由後續 provider 完成。

場資料按 owned rows 保存；resolved boundary arrays 延續既有 runtime 的 replicated
儲存。本批未移除此副本，也未設計其大型磁碟分片。HPC-06A 的記憶體稽核與後續
codec 必須處理這個成本，不能把本批當作完成所有分散式 checkpoint I/O。

## 驗收與結果

[測試程式](../../solvers/cpu/tests/test_body_fitted_accepted_checkpoint.cpp) 重用既有
bifurcation 的 source fixture builder：一個 cubic element、64 nodes，三個 flow ports；
另加入 tracer 的 time derivative、volume source 與 diffusion。這是小型正確性 fixture，
不宣稱 scaling。每組使用穩態、Backward Euler、Backward Euler＋RC outlet 三種模式。
保存第 3 步的 snapshot，fresh runtimes 恢復後完成第 4–6 步。

| 檢查 | 結果 |
|---|---|
| accepted state | 三種模式 × 3／1／2-rank groups，共九種組合；剛恢復及後續三步的 owned flow／scalar fields、boundaries／tractions、outlets、clock、linear iterations、scalar steps 精確一致；共 27 個續跑步驟 |
| 物種質量 | 每個續跑步驟用原生 `TotalMass()` 做 collective 積分，恢復與參考結果精確一致 |
| 局部錯誤 | 只破壞最後一個 local rank 的 identity、clock／dt、ownership、NaN、field 長度、boundary mask／count、traction label 或 iteration count；全群在 validation 階段一致拒絕，fresh fields／counts 保持原值 |
| 群組不同值 | 合法但不同的 replicated boundary value 或 transport step count，在 agreement 階段一致拒絕；沒有先寫進 candidate state |
| lifecycle | initial／active／prepared／closed capture 拒絕；重複 restore 拒絕；prepared abort 後完整 accepted state 不變；transport counter overflow 拒絕 |
| 既有 trial regression | 三 rank job、WORLD 與 1／2 groups 通過原有局部 failure／rollback／retry／COMM_SELF 比較；1／2／3-rank groups 分別報告 33／39／39 項 |
| 既有 staged regression | 三 rank job、WORLD 與 1／2 groups 通過 staged hydraulic／transport 與物種帳目；1／2／3-rank groups 分別報告 45／52／52 項 |
| 既有 VCA runtime unit | 重建後的 `vca_3d_runtime_test` 退出碼 0 |
| 正式 native graph | 0D–3D–0D fixture 執行兩步，對照先前封存的不中斷結果；7 份 CSV、28 rows、204 numeric values，最大絕對差 0 |

新 typed-state 比較門檻事先設為精確相等。既有 failure regressions 沿用原有
COMM_SELF relative 1e-6／zero-reference absolute 1e-12 gate；native graph 的 gate
為 `abs(error) <= 1e-12 + 1e-6*abs(reference)`，非數值欄位與欄列結構須相同。
native 比較驗證 accepted histories，沒有新增完整 3D field export 的比較或宣稱新
CLI checkpoint 已放行。

ASan／UBSan 另以同一新測試 source 建置及三 rank 執行；系統 MPI／PETSc 未經
sanitizer 建置，設定 `ASAN_OPTIONS=detect_leaks=0`，不宣稱 third-party leak audit。
最終作業三個 rank 退出碼均為 0，stderr 全部為空；九組 state／續算比較通過。

## 環境、時間與命令

Linux x86_64 工作站、GCC 11.4、OpenMPI 4.1.2、PETSc 3.15.5 real64/int32。
本機 3 ranks、1 thread、preonly／LU／MUMPS；無固定 affinity。C++17、
`-Wall -Wextra -Wpedantic` 建置，無新警告。restricted build 的 `opal_ifinit errno=1`
是本機 socket 權限訊息；實際 MPI 以所需執行權限完成。

最終測試逐 rank 的 exclusive phase 時間如下，避免把 nested totals 重複相加：

| World rank | assembly s | solver setup s | linear solve s | 明確量測的 communication s | peak RSS bytes |
|---|---:|---:|---:|---:|---:|
| 0 | 7.39510020 | 0.51445183 | 0.03085015 | 0.09098117 | 46,772,224 |
| 1 | 7.27860452 | 0.56008056 | 0.03854161 | 0.33479949 | 45,211,648 |
| 2 | 7.26165187 | 0.56057562 | 0.03858628 | 0.53330921 | 48,525,312 |

每 rank wall 最大 8.83048 s。flow＋transport 一次 typed restore 的計時範圍為
0.000017494–0.000191819 s，包含 local validation、群組 agreement、staging、snapshot
與 publication，不含磁碟 I/O。每 rank 都保留各次 restore 時間。測試同時包含 WORLD
與 split groups，且與其他正確性驗證有執行重疊；這些數字不是獨立 performance
repetitions，不能推論加速或 scaling。磁碟 checkpoint／CUDA allocation 為 N/A。

既有 trial、staged、native jobs 的最大 wall 分別為 3.22067／1.36664／0.41460 s；
逐 rank RSS 與 profile 保存在 `measurements.json`。native output 與失敗測試的
wall 均含初始化與測試工作，不拿它們當單獨 assembly／solve 時間。

```bash
make -C solvers/cpu body_fitted_accepted_checkpoint_test vca_3d_runtime_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
make -C solvers/coupling three_d_trial_failure_test three_d_staged_failure_test \
  iga_multidomain_flow PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
mkdir -p outputs/hpc05/body-fitted-state/final-mpi
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 IGA_PROFILE=1 \
PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps' \
mpiexec --oversubscribe -np 3 python3 scripts/hpc_rank_run.py \
  --output-dir outputs/hpc05/body-fitted-state/final-mpi --expected-ranks 3 --timeout 180 \
  -- solvers/cpu/body_fitted_accepted_checkpoint_test \
  outputs/hpc05/body-fitted-state/final-fixtures
```

每次執行使用新的 directory。sanitizer flags 為 `-O1 -g -fsanitize=address,undefined
-fno-omit-frame-pointer`，執行另設 `UBSAN_OPTIONS=halt_on_error=1`。MPI、unit 與
sanitizer 都使用本機所需執行權限，沒有 cluster login node 大型測試。

證據根目錄：`outputs/hpc05/body-fitted-state/`。包含 `final-mpi`、`sanitized-mpi`、
`trial-mpi`、`staged-mpi`、`native-mpi` 各 rank stdout／stderr／resource／run.json，
VCA unit logs、`native-comparison.json`、`measurements.json`、來源／fixture／binary
雜湊及封存 executable。`acceptance.json` 核對 158 份來源與 182 份證據，
commit 後另以 `commit.json` 核對提交內容。`attempt-1` 是上述 dt 區分修正前的失敗，不列入通過結果。

## 剩餘工作與下一步

下一步把此 typed state 接上有 shape／clock／identity gate 的 3D streaming codec，
加入 bundle 分片與新 MPI 作業載入的驗收。然後完成 1D／3D graph provider wiring、
完整 accepted history prefix、native CLI 的 capture／load／candidate publication、
全作業一致失敗與 accepted clock／pressure guess／donor 恢復，再驗證三個中斷點。

本批未完成 HPC-05C，也未處理 immersed／moving／FSI／不同 rank 數恢復。整份
TODO 仍為 14／38 已完成；05C 保持部分完成，不把 runtime 元件測試當作完整 graph
checkpoint 驗收。
