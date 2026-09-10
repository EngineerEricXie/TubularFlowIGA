# HPC-04B 小切割 gauge 案例

- 狀態：三候選評估完成（LU 通過、兩個 block 候選不收斂）；HPC-04 整體未完成。
- 日期：2026-09-09（美東）；基準 `73889c4`。只新增測試幾何編譯分支與 target，
  production discretization、quadrature defaults、solver options 與資料格式不變。

`immersed_small_cut_solver_test` 重用 static MPI test 全部驗收，將立方體從
[0.1,0.9]³ 改為 [0.3,0.7]³，仍用 3×3×3 background cells；名義 octree depth
由 2 改為 3。角落 cell 的解析 fluid volume fraction 從 0.343 降到 0.001，
本次 quadrature 最小估計值為 0.0011009932472103598。測試斷言存在正值且小於
0.01 的 usable fraction。估計值不是解析幾何體積，不宣稱本次已完成 quadrature
volume convergence；此處評估相同離散算子下的 solver 條件敏感性。

封閉 body-force／wall 案例有 864 physical DOFs 加一個 gauge（row 864）；
Nitsche wall、ghost penalty、原 serial reference、rollback／prepare／commit、
true linear residual、分欄位與守恆 gate 全部保留。

| 候選 | 結果 | 場相對 L2 | peak rank RSS bytes |
|---|---|---:|---:|
| FGMRES／LU／MUMPS | 通過 | 2.83951e-9 | 98,779,136 |
| FGMRES／block-Jacobi／ILU | 原上限不收斂 | 無 accepted field | 96,944,128 |
| FGMRES／block-Jacobi／LU／MUMPS，shift=1e-12 | 原上限不收斂 | 無 accepted field | 99,090,432 |

兩個失敗候選記錄 `DIVERGED_ITS iterations 2000`、reason=-3，所有 ranks 退出 1，
沒有成功的 accepted-state 測試紀錄。沒有增加上限或放寬 tolerance。原較大 fluid
fraction 的 closed／flow／pressure 九候選均通過，故本次小切割失敗是必須保留的
條件敏感性證據，不等同所有 block 預條件器都不可用。

LU 的 pressure relative L2=2.8394485776888653e-9、velocity absolute L2=
3.9258424195645828e-11、scalar absolute L2=1.5371505902501382e-15。原 closed
近零速度／scalar 使用 absolute gate≤1e-10；沒有以接近零參考的相對誤差作判準。
4 次 nonlinear iterations、5 次 KSP iterations，最終 residual=2.7073682537633245e-18；
volume divergence=1.4707832778482977e-19，surface flux=1.4885610850493027e-19。
Rank 0 assembly=41.776107894 s、linear=0.389050169 s；原未啟用小切割編譯分支的
empty-work 2-rank 回歸也通過，全場相對 L2=4.71934e-15。
`make mesh-test` 退出 0、`mesh_core_test: PASS`，紀錄在 `mesh-test.log`。

## 證據與重現

工作站 Intel Core i9-14900KF（可見 16 logical CPUs）、GCC 11.4、Open MPI 4.1.2、
PETSc 3.15.5 real64/int32，2 ranks、OMP／OpenBLAS=1，CPU test 未編入 OpenMP。
同時有背景回歸，時間不作性能排名；RSS 包含 serial reference 與測試成本，
未涉及 CUDA allocation。

本機 `outputs/hpc04/immersed-candidates/small-cut-closed-v1/acceptance.json` 保存
3 作業、6 份 rank report，status=`completed_with_failures`、harness 退出 1。
`small-cut-acceptance.json`／`audit_small_cut.py` 核對一成功、兩失敗、cut fraction
與另兩份 default regression rank report。`small-cut-source.json` 及
`small-cut-v1-binary` 保存來源與 binary 身分。原 regular suite 的 frozen
`static-v1-binary`／`3e3e9b6` 身分保留，沒有用新建置替換其驗收歸因。

```bash
make -C solvers/cpu immersed_small_cut_solver_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
python3 scripts/hpc_immersed_preconditioners.py \
  --binary solvers/cpu/immersed_small_cut_solver_test \
  --output-dir NEW_OUTPUT --modes closed --timeout 1200
```

下一步是小切割 port 條件、更大網格、rank 與其他預條件器評估，以及無背景干擾／
具排程資源的效能驗收；保留目前已驗證 production defaults。
