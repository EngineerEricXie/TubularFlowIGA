# HPC-04A：貼體 graph 的獨立求解器選項

本報告保留各開發階段當時的狀態；目前 HPC-04A 已完成，整合來源與驗收範圍見
[完成稽核](HPC_04A_COMPLETION_AUDIT.md)。HPC-04B／C 仍進行中。

日期：2026-09-09。基準 revision `fa89f05b1e9597d98f148e7d2cf65c2ebbcdb3f0`
加本批修改。**部分完成**：正式 multidomain／bifurcation runner 的貼體 flow、transport
可獨立配置並記錄實際求解器；1D implicit／SNES、其他 standalone、immersed／moving／FSI
及完整巢狀診斷與後端矩陣仍待完成。整份 TODO 保留 38 項目標、14 項已完成。

## 實作與使用

操作範例與優先順序見 [SOLVER_OPTIONS.md](../SOLVER_OPTIONS.md)。
`PetscSolverOptions.hpp` 建立每個 runtime 的精確 options snapshot，依 runtime defaults、
未加前綴的共同選項、domain／role 選項依序覆蓋。prefix 與完整有效選項在 runtime
communicator 內一致後才建立求解器。保留 flag、空白、引號、換行與 option-like value，
不使用顯示字串重新解析；使用旗標傳回來源，來源值保持不變。

KSP 與 PC 先掛入相同 database、設定 prefix，再建立 fieldsplit children。PETSc 3.15
部分子 KSP 只繼承 prefix，因此 SetFromOptions／SetUp／Solve 期間另以 private database
作暫時預設。scope 使用 returning error handler；成功、PETSc returned error、C++ 例外
及拒絕重入皆還原 options／handler stack，再進行跨 rank 錯誤協調。這不能保證任意
第三方 PETSc collective 在任一 rank 內部故障後都能返回；此批證據涵蓋配置錯誤與既有
runtime 故障注入，並未把 HPC-01C 全項勾選完成。

正式 graph 每個 accepted step 寫 `solver_configuration` JSON：domain、role、prefix、
KSP、PC、factor backend、tolerances、iteration limit、最後一次線性 iterations／reason。
最後一次 solve 可因 warm start 已收斂而為 0 iterations；不當作所有 Newton/trial 的
累計值。舊 embedding caller 的空 prefix 保留原介面與 scalable flow 的全域 defaults
行為。原場／history／checkpoint payload 格式不變；來源指紋改變，舊建置 checkpoint
按既有規則拒絕混用。

最終 native source identity：
`0367112220bc6829805f055311372479f0c4295ee2570f3c6ce09d1b5c313880`。

## 環境與重現

本機 `TsungYehLab`，GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real64／int32，MUMPS。
`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`；rank affinity 可用 CPUs 0–15，無固定
逐核 binding。小 fixture 的並行驗收不代表排程器上的 scaling benchmark。

```bash
export PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
make -C solvers/cpu petsc_solver_options_test PETSC_DIR="$PETSC_DIR"
make -C solvers/coupling -j2 iga_multidomain_flow native_graph_checkpoint_fault_test \
  iga_1d_3d_bifurcation iga_1d_3d_explicit runtime_construction_failure_test \
  multidomain_subcommunicator_test PETSC_DIR="$PETSC_DIR"
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 mpiexec -np 3 solvers/cpu/petsc_solver_options_test
python3 scripts/hpc_solver_prefixes.py \
  --binary outputs/hpc04/options/native-handler \
  --reference outputs/hpc05/native-graph/binaries/native-v6 \
  --output-root NEW_EVIDENCE_ROOT --ranks 3
python3 scripts/hpc_native_graph_checkpoint.py \
  --binary outputs/hpc04/options/native-handler \
  --fault-binary outputs/hpc04/options/fault-handler \
  --output-root NEW_RESTART_EVIDENCE_ROOT --ranks 3
```

原生比較腳本自行設定共同 `preonly / lu / mumps`，domain override 保留同一 LU backend，
分別改用 FGMRES／GMRES。reference 為基準 revision 的封存 executable；需保留此 binary
才能重現跨 revision 的精確比較。正向作業以 `hpc_rank_run.py` 保留每 rank exit、RSS、
stdout／stderr；預期失敗使用有 timeout 的直接 MPI 作業，檢查指定退出碼。

## 門檻與結果

未指定 override 時，0D–3D、1D–3D、兩物種及雙 3D island 的所有輸出檔與各步 field
payload bytes 必須精確相同。指定 override 時，場及有量綱 history 同時滿足 L2 與逐項
`|difference| <= 1e-12 + 1e-6 * |reference|`，數值有限；IDs、step、iterations 精確。
原有 Newton、coupling 與 species amount gates 不變。已正規化的 species residual 是
無量綱診斷，另按 case 原始 `relative=1e-6`、`absolute=1e-10`、`reference_amount=1e-5`
檢查原始 gate 與絕對診斷差，不再對接近 roundoff 的正規化 residual 做第二次相對誤差。

| 驗收 | 結果與證據（`outputs/hpc04/options/` 下） |
|---|---|
| 2／3／4 ranks 原生 prefix | `handler-ranks2`、`handler-suite`、`handler-ranks4` 各 14 作業通過，共 36 正向、6 預期 exit 1；預設結果 bytes 完全相同 |
| 不同 domain 與 role | 雙 island 分別 FGMRES／GMRES，species 的 flow FGMRES／transport GMRES；實際 KSP／PC／MUMPS 與預期相同，reason 全為正 |
| 最大全場相對 L2 | 2 ranks：`1.953e-12`；3 ranks：`1.587e-12`；4 ranks：`1.390e-12` |
| 不支援選項 | 多 rank scoped LU backend `petsc`、不存在的 scoped transport KSP 均正常 exit 1；其後健康新作業結果與 reference 精確相同 |
| 單 rank 原生 prefix | `handler-rank1` 的 flow／species／雙 island 共 9 作業通過；預設結果精確一致，override 場相對 L2 最大 `4.732e-13` |
| 最終 checkpoint 回歸 | `checkpoint-handler/acceptance.json`：63 作業通過，41 正向、22 預期失敗；full／save／resume、故障世代恢復、身分拒絕、signal stop 與恢復保持既有驗收 |
| 單元 | `unit-handler-mpi`：world 3 ranks 及 split 1／2 ranks；block-Jacobi、additive fieldsplit、Schur FULL/A11 的子 solver override 與 MUMPS backend 正確，精確 SPD 解誤差低於 `1e-10` |
| 建構／清理及子群組 | `construction-handler`、`subgroups-handler-flow`、`subgroups-handler-species` 通過；world 3 及 split 1／2 ranks 故障注入後釋放 retained references、重試場一致、scalable defaults rollback 與使用者 override 保留 |
| ASan／UBSan | `handler-sanitized-suite` 的 14 作業通過（12 正向、2 預期 exit 1），default bytes 一致、override 場與 history 通過相同門檻，無 sanitizer 診斷 |
| scope 邊界 | 同一測試覆蓋單 rank prefix／options 分歧、非法與超長 prefix、長來源 key、實際未知 KSP、單 rank returned error／C++ 例外、重入拒絕與健康重試；caller handler 計數探針驗證成功／失敗／例外後還原 |

最終集合共 132 個 MPI 作業（102 正向、30 預期失敗），正向 288 份 rank reports。
`outputs/hpc04/options/acceptance.json` 記錄 source／evidence SHA256、所有作業及報告；
同目錄 `audit_final.py` 重新核對退出碼、stdout／stderr hashes、RSS 與 sanitizer field
catalog。較早的 `checkpoint-regression` 等通過紀錄不混入最終集合。
ASan／UBSan 以 `-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer` 編譯，
`ASAN_OPTIONS=detect_leaks=0:halt_on_error=1`、`UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1`；
關閉第三方函式庫的 leak scanning，因此不聲稱完整第三方 leak audit。完整編譯 argv
保存於 `sanitized-handler-command.json`。所有 affected native／standalone 目標重建，
沒有新增 compiler warning；Make 包裝器在 sandbox 印出的 Open MPI socket probe 訊息
不影響編譯，實際 MPI 作業均在允許 MPI 的環境執行。

三 ranks 的 override 作業下表為每 rank 已計時 phase 的最大 inclusive 秒數與最高 RSS。
通信／輸出可能巢狀，不可相加；未包入 phase 的工作仍是 unscoped。資料來自
`handler-metrics.json` 與原始 profile。這是小案例成本紀錄，未作預條件器效能優劣結論；
本批未執行 CUDA，CUDA peak allocation 不適用。

| 案例 | assembly s | solver setup s | linear solve s | communication s | output s | peak rank RSS MiB |
|---|---:|---:|---:|---:|---:|---:|
| flow | 0.3262 | 0.0210 | 0.0017 | 0.4884 | 0.4441 | 43.92 |
| species | 0.3510 | 0.0324 | 0.0038 | 0.6597 | 0.5601 | 45.25 |
| 雙 island | 0.7332 | 0.0366 | 0.0078 | 0.9559 | 0.7589 | 46.12 |

## 開發中發現及剩餘工作

早期 `native-suite` 將 warm-start 的 0 iterations 當失敗，`native-suite-v2` 對近零的
normalized residual 重複取相對差；修正驗收分類後，物理場與原守恆 gates 都通過。
`native-final-suite` 揭露實際錯誤：預設 PETSc handler 在未知 transport KSP 時令其他
ranks SIGABRT（exit 134）。加入 scoped returning handler 後，最終 2／3／4 ranks 均
exit 1，未放寬預期退出碼。這些失敗日誌仍保留，不列入最終通過數。

下一步接入 1D implicit／SNES 與其 MUMPS 要求，再處理其餘 runtime 的 prefix 和完整
nested solver 診斷。HPC-04B/C 尚須用固定數值門檻驗證不同網格與後端的成本、失敗及
負收益組合。大型／跨節點排程驗收仍需可用 cluster allocation；目前沒有此項通過證據。
