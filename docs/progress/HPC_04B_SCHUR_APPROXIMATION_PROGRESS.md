# HPC-04B Schur 近似評估

- 狀態：16-element 小案例評估完成，包含不收斂候選；HPC-04B 整體仍進行中。
- 日期：2026-09-09（美東），基準 `42d4d39`。Native frozen binary 仍為
  `01aa4c8`；未更動 production solver 或物理／資料介面。
- 環境：工作站回報 Intel Core i9-14900KF，執行環境可見 16 logical CPUs；
  PETSc 3.15.5 real64/int32、Open MPI 4.1.2、GCC 11.4，2 ranks，OMP／OpenBLAS=1。
  同時有其他功能回歸，且本 suite 開啟 KSP view，不作無干擾效能排名。

沿用 [方管候選](HPC_04B_DUCT_CANDIDATES_PROGRESS.md) 的 16-element、175-node、
兩步強耦合案例及全部原 Phase 9 validator／checkpoint 場門檻。
先前 A11 Schur／LU 累計 62,351 次線性迭代，故加入 `selfp`：PETSc 實際 view
確認使用以 A00 對角倒數形成的 assembled Schur approximation。只改候選選項，
保留 FGMRES、rtol=1e-12、maximum iterations=5000。

| 候選 | 速度 PC | 壓力 PC | 結果 | 累計線性迭代 | 最大場相對 L2 |
|---|---|---|---|---:|---:|
| lu | 全系統 LU／MUMPS | 同左 | 通過 | 42 | 參考 |
| selfp-lu | LU／MUMPS | LU／MUMPS | 通過 | 1,967 | 7.08198e-12 |
| selfp-gamg | GAMG | Jacobi | 第一步不收斂 | 5,000（失敗 solve） | 無 accepted field |
| selfp-gamg-lu | GAMG | LU／MUMPS | 通過 | 12,484 | 3.20496e-11 |

`selfp-gamg` 於第一個 nonlinear iteration 達 KSP 上限，reason=-3，reported residual
3.3675428753162122e-10；所有 ranks 退出 1，未發布 accepted output。
沒有提高迭代上限或放寬 tolerance 使它通過。改用壓力 LU 後可以通過，是這個具體
組合的證據，不推論所有壓力 Jacobi 或 GAMG 都不可用。

三個成功作業各完成 42 次線性 solve、兩個 accepted steps，通過原 pressure／flow／
3D mass／0D constitutive／global volume gates。實際 child PC view 與配置選項相符。
相較先前 A11／LU，selfp／LU 在相同小案例減少迭代，但仍不能據此推論更大網格效能。

## 資源與證據

下列為各 rank 最大 exclusive phase seconds，並非同一 rank 的可相加時間：

| 候選 | assembly | setup | linear solve | communication | output | peak rank RSS bytes |
|---|---:|---:|---:|---:|---:|---:|
| lu | 34.1998 | 1.05455 | 0.0526398 | 1.11290 | 0.165824 | 60,354,560 |
| selfp-lu | 34.0768 | 0.303926 | 2.06112 | 1.06444 | 0.166389 | 63,365,120 |
| selfp-gamg-lu | 34.4351 | 0.318276 | 8.99566 | 1.04310 | 0.175821 | 55,906,304 |

不收斂的 selfp-gamg peak rank RSS 為 51,429,376 bytes；失敗作業不列 accepted-step
效能比較。未涉及 CUDA allocation。

`outputs/hpc04/pde-candidates/selfp-small-v2/acceptance.json` 保存 4 個作業與 8 份
rank report，包括失敗的 resource、stderr tail；`audit.json`／`audit_selfp.py` 檢查
三個通過、一個失敗，核對實際 Schur 近似、兩個 child PC、迭代上限、無失敗輸出，
並保存 acceptance hash。`selfp-small-v1` 是新增失敗診斷保存前的三候選評估，
不與 v2 相加當作額外完成範圍。

Harness 新增 `--candidates` 與 `--solver-view`；LU 參考始終先跑，沒有指定候選時
保留原四候選預設。候選不收斂時，整體 status 是 `completed_with_failures` 且退出 1，
不把「評估已結束」寫成數值全通過。

```bash
python3 scripts/hpc_duct_solver_candidates.py \
  --binary outputs/hpc04/immersed/native-v4-binary \
  --output-dir NEW_OUTPUT \
  --candidates selfp-lu selfp-gamg selfp-gamg-lu --solver-view
```

## 剩餘工作

原四候選的 128-element、539-node、4-rank sweep 已啟動，證據目錄
`outputs/hpc04/pde-candidates/refined-four-v1`，每個候選 timeout=1200 秒，尚未完成。
該作業使用加入失敗 stderr 保存前的 interim harness，report 記錄其獨立 harness hash；
不能歸為本次最終 harness revision 的實測。
還須不同 rank 的同網格場比較、更多 pressure reference／nullspace 與 immersed 條件，
以及排程資源的大案例與無干擾 timing。完整 moving 回歸已在 rigid wall trace gate 失敗，定位進度見 [immersed options 回歸](HPC_04A_IMMERSED_OPTIONS_PROGRESS.md)。
