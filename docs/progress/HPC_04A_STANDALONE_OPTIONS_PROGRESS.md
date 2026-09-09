# HPC-04A：Body-fitted standalone／VCA 的 solver prefix

日期：2026-09-09。基準 commit `6393c7dcbf11833a497458156920405a2a3f6166` 加本批修改。
**部分完成**：`iga_navier_stokes`、`iga_solve` 及 VCA transport 已接入獨立配置；
immersed／moving／FSI、完整巢狀診斷及後端矩陣仍待完成。整份 TODO 仍有 38 項，
14 項已完成。操作見 [SOLVER_OPTIONS.md](../SOLVER_OPTIONS.md)。

## 變更

Standalone flow 以 Navier–Stokes system 名稱產生 `..._flow_` prefix；沒有 system 的
legacy boundary case 使用 `domain_flow_flow_`。VCA 的 transport 以自身 system 名稱
產生 `..._transport_`。原 `TransientFlowRuntime`／`TransientTransportRuntime` 的
private snapshot、既有 typed defaults、後端能力檢查及錯誤協調沿用上一批實作。

一般 transport CLI 的 `TransportPetscObjects` 持有 options owner，活到 KSP、vectors
及 matrices 清理完成。SetFromOptions、SetUp、SetUpOnBlocks、Solve 都使用 private
scope，使用旗標回傳來源；既有 GMRES／block-Jacobi 預設與 warm start 不變。
每一步記錄有效 prefix、KSP、PC、factor backend、iterations／reason，stdout 寫入
仍位於 collective error boundary。這是最後一次線性呼叫的值，不是 Newton 總和。

兩個 CLI 原先會把單 `-` 開頭的 PETSc options 當成 legacy positional argument；
現先略過 PETSc key 及 optional value，由 PETSc 在初始化時解析。`--` application
options、既有 positional 介面與 application option consumption 保留；位置參數中的
負數不會被誤當作字母開頭的 PETSc key。具值 option 後面的負數亦可交給 PETSc。

沒有變更 `.ntiga`、場、VTK、checkpoint 或 VCA 帳目格式；僅 CLI source／Makefile
改變，共用 runtime header 與 native graph source identity 未變。

## 驗收及限制

本機 `TsungYehLab`，GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real64／int32，MUMPS。
本次 CPU CLI 為 Makefile 預設的非 OpenMP 建置；runtime 明確報告
`openmp_compiled=0`。OpenBLAS／OMP 環境各 1 thread，允許 CPU 0–15，沒有固定逐核
binding。建置無新增 compiler warnings。本批未作 CUDA、大型或跨節點驗收。

`outputs/hpc04/standalone/final-suite-v3/acceptance.json` 的 **27 個 MPI 作業全部符合
預期：18 正向、9 預期 exit 1，27 份正向 rank reports**。使用既有 healthy VCA flow、
prescribed transport、velocity-series transport，各 1／2 ranks。執行前核對來源 manifest
記錄的全部 input SHA256；immutable input 不在此批修改。reference binary 在修改前
封存，SHA256 保存於 acceptance。

| 檢查 | 結果 |
|---|---|
| 未指定 domain override | 所有場、pressure、VTK／index、checkpoint、VCA metadata／manifest bytes 與 reference 精確一致 |
| 前綴 FGMRES override | 有效 prefix／FGMRES／LU／MUMPS 記錄符合預期；40 個 text、pressure、完整 flow／transport checkpoint 場比較通過 |
| 場數值門檻 | L2 difference ≤ `1e-12 + 1e-6 × reference L2`，相同長度且全為有限值；最大相對 L2 `2.766e-12` |
| VCA／checkpoint JSON | 14 個 metadata 比較通過；keys、型別、IDs、strings 精確，浮點逐值門檻 `1e-12 + 1e-6 × |reference|` |
| 未知前綴 KSP | 六個案例均在預期的 flow solver options／transport KSPSetFromOptions 邊界正常 exit 1 |
| 不支援的多 rank LU | 三個 2-rank 案例明確拒絕 scoped `pc_factor_mat_solver_type=petsc`，診斷為 factor backend availability |
| VCA smoke | `vca-override.log`：7 個 solver 子作業（含 3 個 2-rank），flow FGMRES／oxygen transport GMRES；原質量帳目、RCR、oxygenator source、checkpoint／restart assertions 全部通過 |

VCA smoke 原有門檻直接保留，例如 oxygenator mass increase 與 source 差不超過
`2e-8`。未將 Krylov iteration 變少當作效能收益。下列是每個 override 作業的最大
rank inclusive phase 秒數與最高 rank RSS；巢狀通信／輸出不可相加，其他未計時工作
仍屬 unscoped。VCA smoke harness 未另外量測整體 peak RSS；不以它作 memory benchmark。

| 案例 | ranks | assembly s | setup s | linear solve s | communication s | output s | RSS MiB |
|---|---:|---:|---:|---:|---:|---:|---:|
| VCA flow | 1 | 0.09635 | 0.00740 | 0.00031 | 0.00028 | 0.00769 | 43.02 |
| VCA flow | 2 | 0.10758 | 0.00750 | 0.00046 | 0.03717 | 0.01074 | 43.60 |
| prescribed transport | 1 | 0.00218 | 0.00071 | 0.00018 | 0.00016 | 0.00468 | 39.37 |
| prescribed transport | 2 | 0.00228 | 0.00094 | 0.00053 | 0.00400 | 0.00530 | 39.48 |
| series transport | 1 | 0.00428 | 0.00091 | 0.00016 | 0.00024 | 0.00514 | 39.29 |
| series transport | 2 | 0.00508 | 0.00129 | 0.00038 | 0.00404 | 0.00484 | 39.66 |

## 重現

```bash
export PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
make -C solvers/cpu -j2 iga_solve iga_navier_stokes PETSC_DIR="$PETSC_DIR"
python3 scripts/hpc_standalone_solver_prefixes.py \
  --output-root NEW_ROOT --reference-dir outputs/hpc04/standalone \
  --flow-manifest outputs/hpc01/flow-step/cli-final/summary.json \
  --transport-manifest outputs/hpc01/transport-cli/cli-final/summary.json
```

此腳本重播既有 healthy acceptance corpus；需保留 manifest 中的輸入檔與兩個封存
reference binary。也可先用 `hpc_flow_step_regression.py`／`hpc_transport_cli_regression.py`
重新產生相同 corpus，再傳入新 manifest。完整 evidence／source hashes 與核對方式
位於 `outputs/hpc04/standalone/acceptance.json` 及 `audit_final.py`。

早期 probe 用 PREONLY 當一般 transport 基線，遇到既有 nonzero initial guess 限制
而失敗；改回來源 manifest 原已驗證的 GMRES baseline，程式沒有改 warm-start 策略。
另一個 probe 假定 `rtol=-1` 必須被 PETSc 拒絕，但本機 build 接受它；因此最終錯誤
驗收使用已確認不支援的 KSP type／factor backend，不聲稱 CLI 新增數值參數驗證。
這些早期失敗及假設修正均保留日誌，不列入最終通過數。本批未另跑 sanitizer，沒有
新增 sanitizer／leak audit 結論。

下一步接入 immersed static／distributed／transient 及 moving／FSI 的 solver options，
補齊完整 nested diagnostics 和後端矩陣。HPC-04A 保留部分完成；HPC-04B/C 的候選
預條件器及中／大案例數值與成本驗收亦未被本批小案例取代。
