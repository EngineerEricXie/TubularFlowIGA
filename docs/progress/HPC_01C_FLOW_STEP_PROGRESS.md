# HPC-01C：CPU flow 每步輸入與 VCA 耦合錯誤

日期：2026-09-08。狀態：本批驗收通過；HPC-01C 與整份清單仍未完成。

## 實作

原生 `iga_navier_stokes` 在進入求解或下一個 collective 前，仍有未協調的
本地配置、waveform 讀取與 VCA 結果處理。本批補上以下邊界：

- VCA face catalogue、pressure traction catalogue 的配置與建立。
- VCA flow／transport required-node layout 檢查、checkpoint identity 與
  history writer 建構，以及 history 輸出的 root 錯誤。
- 初始 VCA 邊界與每步配置、waveform、幫浦入口及物種邊界準備。
- 原生 CLI 呼叫本地 `SetTrialBoundaryConfiguration` 的錯誤。
- VCA port 結果配置、物種 flux／concentration 複製、transport budget 的
  本地計算與 previous-mass 更新，以及 circuit advance／history bookkeeping。

MPI／PETSc collective 保持在 local-stage callback 之外。新增的準備函式
以 optional 在 callback 內建構結果，並用 compile-time assertion 確認返回
時的 move constructor 不會拋例外，避免在共同驗證後再新增配置失敗點。

另將條件運算式中的臨時空物種清單改成既有空清單參照：兩個分支均為 lvalue，
不再因一個分支是臨時 vector 而在進入 port collective 前複製所有物種名稱。
velocity gather 也改成明確的前一步操作，再將結果參照傳入 transport。

維持原有時間順序：暫態 waveform 使用目前物理時間；VCA inlet 使用前一步
時間；初始 VCA inlet 使用 0。steady 路徑保留配置複製行為，legacy 路徑保留
原有預設。數值算式、求解器、時間步順序、接受條件與檔案格式均未變更。

## 驗收

| 驗收 | 結果與範圍 |
|---|---|
| 原生函式的 MPI 故障測試 | world 3 ranks 與獨立 1+2 groups，各 22 個故障及健康重試，共 66／66 通過 |
| 新原生 VCA CLI regression | 13 筆觀察、20 份 rank reports；7 筆健康／刻意停止、6 筆預期失敗 |
| 既有 flow input regression | 51 筆觀察、94 份 rank reports、32 份指定場比較通過；保留 legacy 跨 rank 壓力的 4 筆未通過觀察 |
| 既有 flow output regression | 46 筆觀察、83 份 rank reports、68 份指定場比較通過，含 OpenMP 路徑及舊版錯誤摘要觀察 |
| 同版 writer 故障回歸 | 33 個故障與 33 次重試通過 |
| 同版 VCA smoke | flow／transport／reservoir、守恆與單／雙 rank checkpoint/restart gates 通過 |

66 個函式故障測試直接使用原生 CLI 的實作，涵蓋最後一個 rank 的 temporal
function／table 錯誤、缺失 flow system、無效 reference flow、物種映射／邊界
錯誤、非有限 inlet time、缺失 port 結果、無效 dt、空／重複 outlets、非有限
流量、無效 aggregation epsilon，以及 reservoir volume／species mass 失敗。
每次均檢查共同診斷與後續健康執行。三種 group 使用不同濃度，單 rank group
另做額外準備步驟，使工作次數不相同，驗證 communicator 獨立性。

健康檢查確認 waveform 在 `t=0.01` 產生預期 scale、VCA inlet 保持前一步時間、
物種輸入及 reference-profile scale 正確，並檢查閉迴路 volume／species
守恆。Circuit 失敗後使用重新建構的 circuit 與 history 做健康重試；這不是
整個 VCA macro-step 的原地 rollback 證據。

原生 CLI 使用既有 schema 允許的負 infusion rate，分別令 reservoir 在第一步
或第二步耗盡物種質量。原有方程拒絕負質量，修正版在
`flow VCA circuit advance` 共同退出 1，沒有最後成功摘要或後續場輸出。
第二步失敗後，flow／VCA 的兩份 metadata 與兩份 state 檔案，均與相同配置
刻意停在第一步的執行逐位元組相同；單／雙 rank 均驗證了這四份檔案。
這證明該故障位置沒有覆蓋最後接受步的 checkpoint，不代表一般中斷下的
多檔原子發布協議已完成。

其他新 CLI 案例包括健康完整執行、不指定場輸出、刻意停止，以及 root history
目標為目錄／FIFO。History 型別錯誤在 `flow VCA history output` 共同退出，
沒有成功摘要；此時模擬步已接受，完整第二步 checkpoint 保留。

新 CLI 的 6 份場比較包含 2 份參考自比及 4 份候選比較，最大 relative L2
為 `1.5869627530656535e-13`。既有 input／output regression 的指定比較最大值
為 `1.8156697926380728e-13`。CPU 門檻仍為 relative `1e-6`，零參考 absolute
`1e-12`；legacy 限制見 [flow 輸入報告](HPC_01C_FLOW_INPUT_PROGRESS.md)。

## 命令、環境與證據

本機 GCC 11.4、OpenMPI 4.1.2、PETSc 3.15.5 real64／Int32、MUMPS；
原生 CLI 使用 PREONLY＋MUMPS LU、`ksp_rtol=1e-12`。一般 OMP／BLAS 固定 1，
既有 output regression 的 OpenMP 案例每 rank 使用 2 threads。
小型正確性作業有部分重疊，沒有新增效能或跨節點擴展性宣告。

```bash
make -C solvers/cpu iga_navier_stokes iga_navier_stokes_openmp \
  flow_step_failure_test flow_output_failure_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  timeout --kill-after=5s 90s mpiexec --map-by core --bind-to core -np 3 \
  solvers/cpu/flow_step_failure_test \
  outputs/hpc01/assets/vca-cli/tubularflowiga-vca-3d-smoke \
  outputs/hpc01/flow-step/unit-final
python3 scripts/hpc_flow_step_regression.py \
  --case-dir outputs/hpc01/assets/vca-cli/tubularflowiga-vca-3d-smoke \
  --reference-binary outputs/hpc01/flow-step/iga_navier_stokes-before \
  --output-dir outputs/hpc01/flow-step/cli-final
```

重跑須使用新目錄。`outputs/hpc01/flow-step/` 保存來源／binary 封存、
build logs、`unit-final.log`、`output-unit.log`、三組 CLI 的 `summary.json`、
`vca-smoke/result.json` 與 `acceptance.json`。較早的 `unit/` 與 `cli/` 也通過
當時範圍，最終版另增加不等工作次數的 communicator 檢查及 checkpoint 比較。

## 限制與接續工作

本批協調可返回的本地錯誤，使仍存活的 ranks 共同退出；沒有增加整體 VCA
步驟的 rollback。Flow／transport 仍按原先順序提交，部分 ranks 的 circuit
可能在共同發現錯誤前已更新，本地 history 配置失敗也可能留下部分記憶體狀態。
CLI 會結束該次作業，不能直接把這個 helper 當作可原地重試的交易介面。

輸入仍須保持不變。沒有逐一注入所有 allocation point；其他 CLI／adapter／
executor 與 runtime 邊界仍須核對。HDF5 close-error、geometry report 與
history writer 的關檔檢查，以及完整 checkpoint 發布協議仍待完成。
不將局部錯誤協調當作 MPI 程序故障恢復或跨節點能力證據。
