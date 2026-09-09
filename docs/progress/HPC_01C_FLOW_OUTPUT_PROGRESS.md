# HPC-01C：CPU flow 輸出與失敗後清理

日期：2026-09-08。狀態：本批驗收通過；HPC-01C 與整份清單仍未完成。

## 修正範圍

原生 `iga_navier_stokes` 的最終輸出原先在成功摘要之後才執行。場檔寫入失敗
時，root 亦可能未歸還 PETSc 借用陣列，gather 與正常清理的返回值未檢查。
PVD／速度序列索引的 root 例外沒有先同步，就可能進入 runtime 解構。

現在 `WriteFlowOutput` 從 state 取得 communicator，先驗證 `4*nodes` 的
大小與整數範圍，再分別協調 gather create／begin／end 的返回值。暫存
scatter／Vec 由既有 `PetscGatherObjects` 管理；root 使用 `PetscReadArray`
歸還借用陣列，包含例外退出路徑。正常 cleanup 的返回值亦需通過。

速度與壓力文字檔明確關閉後才檢查寫入結果。文字場、pressure sidecar、VTU、
PVD、velocity manifest，以及 Bezier geometry report／VTKHDF 初始化目標的
既有非 regular file，在開啟前拒絕。root 的場輸出與初始化錯誤透過共同
local stage 報告；輸出路徑配置、時間索引 bookkeeping 與最後索引寫入也已協調。

`navier_stokes_v2 seconds=...` 移到最終場與索引 writer 成功返回後才列印。
`seconds` 保留原先的量測區間，沒有因摘要移動而納入原先未包含的最終輸出時間。
方程、求解器、場資料格式、輸出頻率與數值容許值均不變。

## 驗收

| 驗收 | 結果 |
|---|---|
| 直接使用原生 CLI writer 的 MPI 測試 | world 3 ranks、獨立 1+2 groups 各 11 個故障及重試，共 33／33 通過 |
| 新原生輸出 CLI regression | 46 筆觀察、83 份 rank reports 通過指定條件；含 12 筆健康、33 筆修正版預期失敗與 1 筆舊版失敗觀察 |
| 場比較 | 68 份，含 8 份參考自比、8 份舊版跨輸出格式比較及 52 份改版候選比較 |
| HDF5 實際讀取 | 兩次比較，每次讀取時間、Points、Connectivity、Offsets、velocity 與 pressure datasets |
| 輸出索引 | 7 組 PVD／VTU 可解析、時間與檔案引用正確；有逐步輸出的 5 組 velocity CSV 時間與引用正確 |
| 原有 flow input regression | 51 筆觀察、94 份 rank reports 與 32 份指定場比較通過；保留 4 筆 legacy 跨 rank 壓力未通過觀察 |
| 原有 VCA smoke | 流場／transport／reservoir、守恆及單／雙 rank checkpoint/restart gates 通過 |

新 CLI 健康案例涵蓋 VTU、VTKHDF、僅最終輸出、不輸出，以及 OpenMP binary
的 1／2 ranks、每 rank 2 threads。場最大 relative L2 為
`1.8156697926380728e-13`，仍使用 CPU relative `1e-6` 與零參考 absolute
`1e-12`。小案例不構成效能或跨節點驗收，OpenMP 比較亦不代表已量測加速。

故障案例涵蓋初始、第一步、第二步與最終速度／壓力／VTU，PVD、velocity
manifest、geometry report 與 VTKHDF 目標為目錄／FIFO，以及實際檔案大小
限制。每個修正版失敗 rank 均核對指定共同診斷、退出 1、無 timeout，且沒有
最後成功摘要。步中輸出失敗保留前一步產物，沒有接續最終輸出。

原生 writer 測試另涵蓋最後一個 rank 的 layout／整數上限錯誤、缺失 mesh、
缺失 HDF writer，以及無 owned rows 的非 root ranks。故障後核對來源 Vec
handle 與值不變，再重新輸出並檢查速度、壓力與 VTU。以各 group 不同的數值
檢查獨立 communicator 沒有混用 world。

修改前 binary 在最終速度輸出是目錄時退出 1，卻已印出成功摘要；這筆實際
觀察保留為 `old-final-failure`。修正版對應案例共同報錯且沒有成功摘要。

## 測試修正與限制注入

首輪 `cli/` 在第一個參考執行後遇到 harness 的 `str`／`Path` 型別錯誤；
沒有當作產品失敗，修正腳本並保留原紀錄與舊腳本。

較早 `cli-verified/` 的雙 rank 全程序 `RLIMIT_FSIZE=64` 在輸出前逾時；
`limit-tcp-probe/` 改用 TCP MPI 傳輸仍逾時。兩次原始退出 124 與日誌保留。
此注入同時限制 MPI 啟動所用檔案，未能證明進入目標 writer；本批未定位
MPI 內部的精確等待點，也沒有宣稱這兩次通過。

多 rank 寫入故障改在直接使用原生 writer 的測試中，完成 MPI 初始化後才對
各 group root 設定 `RLIMIT_FSIZE=1`，並忽略 `SIGXFSZ`。實際文字場截斷於
1 byte，共同錯誤返回後恢復原限制與 signal handler，再驗證來源狀態及重試。
此測試在三種 group 配置都通過。原生 CLI 的單 rank 64-byte 注入也通過；
其 stdout／stderr 經無限制父程序的 pipe 保存，避免錯誤日誌本身被截斷。

最終完整 CLI 證據為 `cli-final/`，沒有使用 `--only` 選項，亦未放寬 timeout
或數值門檻。上述替代測試證明的是 writer 的 rank-local I/O 失敗協調，
不證明 MPI 能在全程序檔案限制下啟動。

## 命令與證據

環境：GCC 11.4、OpenMPI 4.1.2、PETSc 3.15.5 real64／Int32、MUMPS。
原生 CLI 為 PREONLY＋MUMPS LU，`ksp_rtol=1e-12`；BLAS 固定 1 thread。
一般案例 OMP／assembly threads 為 1，OpenMP 案例設為 2。
CLI child timeout 60 秒、rank report rendezvous 15 秒、job timeout 90 秒、
kill grace 5 秒。VCA smoke 使用獨立 TMPDIR、外層 180 秒。

```bash
make -C solvers/cpu iga_navier_stokes iga_navier_stokes_openmp \
  flow_output_failure_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  timeout --kill-after=5s 90s mpiexec --map-by core --bind-to core -np 3 \
  solvers/cpu/flow_output_failure_test outputs/hpc01/flow-output/unit-final
python3 scripts/hpc_flow_output_regression.py \
  --case-dir outputs/hpc01/assets/vca-cli/tubularflowiga-vca-3d-smoke \
  --reference-binary outputs/hpc01/flow-output/iga_navier_stokes-before \
  --output-dir outputs/hpc01/flow-output/cli-final
```

以上目錄已存在；重跑須換新目錄。`outputs/hpc01/flow-output/` 保存
`acceptance.json`、來源／binary 封存、build logs、`unit-final.log`、
`cli-final/summary.json`、`input-regression/summary.json`、
`hdf-validation.json`、`index-validation.json` 與 `vca-smoke/result.json`。
HDF reader 來源與實際命令亦保存在同一證據目錄。

## 尚待完成

失敗可能留下部分輸出，沒有新增原子發布、檔案鎖或持久化保證。既有目標型別
檢查不防止檢查後被並行替換。HDF5 close-error 與 geometry report 關檔檢查、
VCA history writer 的完整寫入檢查仍待補齊；成功摘要目前也不證明所有解構／
MPI finalization 已完成。

flow 的 port catalogue、VCA identity／history 配置、每步本地輸入與耦合
bookkeeping，以及其他 CLI／adapter／executor 邊界繼續由 HPC-01C 追蹤。
root 仍收集完整場；分散式 I/O、完整 checkpoint 發布與跨節點驗收分別屬於
HPC-06、HPC-05、HPC-09，不能由本批小案例推論完成。
