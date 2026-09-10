# HPC-04B immersed 候選驗證

- 狀態：進行中，原切割幾何的 closed／flow／pressure 九候選均通過。
- 日期：2026-09-09（美東），基準 `d3a08e6`。本次只擴充既有 static MPI test
  的 layout 診斷及候選 harness；production solver／operator／資料格式未變。

`scripts/hpc_immersed_preconditioners.py` 使用原
`immersed_distributed_static_flow_test` 的切割立方體（3×3×3 background cells，
surface 範圍 [0.1,0.9]³，depth-2 volume quadrature、Nitsche wall、ghost penalty），
每個 rank 各建立原序列參考，再跑分散式 runtime。此幾何有切割元素，但不是任意
極小 volume-fraction 的 sliver 驗收；仍須另行覆蓋更差的切割條件。

三種 mode 保留原定義：closed 以 body force 與封閉壁形成壓力場、使用 gauge；
flow 具有兩個流量控制 port 及 gauge；pressure 將出口改為壓力條件，無 gauge。
測試新增實際 physical／total DOFs、gauge row 與 port 數的輸出，並斷言 gauge
是否符合 mode。Gauge／flow-controller scalar row 追加在物理 DOFs 之後，因此
不能直接套用方管均勻 block-size=4 的整體 fieldsplit。

候選僅 override `domain_test_flow_`，serial reference 維持原 family options：

| 候選 | 分散式 PC |
|---|---|
| lu | LU／MUMPS |
| block-ilu | block-Jacobi，子 KSP preonly／ILU |
| block-lu-shift | block-Jacobi，子 KSP preonly／LU／MUMPS，nonzero pivot shift=1e-12 |

共同使用 FGMRES、rtol=1e-12；原 test nonlinear absolute=1e-14、relative=1e-11
及既有 gate 保留。Pivot shift 只作用於 PC，不更改 operator 或 true residual。
沒有盲目對包含 gauge／controller 的整個 saddle-point 系統套標量 AMG。

## 已完成證據

`outputs/hpc04/immersed-candidates/closed-v1/closed-lu` 兩個 ranks 均退出 0，
physical DOFs=864、total=865、gauge row=864、無 ports，owned rows=432／433。
原 transaction rollback、failed prepare／abort、commit、true linear residual、
分欄位誤差與守恆 gate 全部通過。

- 全場相對 L2：4.6060734556832725e-14。
- 壓力相對 L2：4.5894657961689117e-14。
- 速度絕對 L2：8.1265829941743591e-15；scalar 絕對 L2：1.0412522250564905e-18。
- Volume divergence：-2.3515603540790086e-17；surface flux：-2.8736721256901933e-17。
- 3 次 nonlinear iterations、3 次 KSP iterations；最終 residual=7.9914149380046672e-16。
- 最大 rank assembly=37.448733319 s、linear=0.370319830 s；peak rank RSS=97,738,752 bytes。

封閉解的速度／gauge multiplier 參考近零，原 gate 使用 absolute L2≤1e-10；
這裡 velocity relative L2 約 36.67，不把它當作有意義的相對精度，也沒有為本次
候選更改 gate。非封閉 mode 仍要求原 velocity／scalar relative L2≤1e-6。

環境為工作站 Intel Core i9-14900KF（可見 16 logical CPUs），PETSc 3.15.5
real64/int32、Open MPI 4.1.2、GCC 11.4，OMP／OpenBLAS=1、CPU test 未啟用 OpenMP。
時間及 RSS 包含既有測試工作，另有背景回歸；不能作候選效能排名。
未涉及 CUDA allocation。Source／binary hashes 在 `source-binary.json`，
`static-v1-binary` 保存本次測試 binary；舊 `01aa4c8` static 結果仍使用當時
`distributed-final-v2/summary.json` 的 binary hash，不能換成新 binary 歸因。

## 重現與剩餘工作

```bash
make -C solvers/cpu immersed_distributed_static_flow_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
python3 scripts/hpc_immersed_preconditioners.py --output-dir NEW_OUTPUT
```

`--modes closed` 或 `--modes flow pressure` 可分批執行；預設每個候選 timeout=600 秒，
並記錄所有 rank 的 resource、log hash、diagnostic tail，失敗不會冒充成功。
LU 參考失敗會停止該次 evaluation；其他候選失敗則保存並繼續，整體返回 1。

`closed-v1` 的三候選、`ports-v1` 的 flow 三候選已全部通過（6 作業、12 rank reports）。
封閉 block-ILU／block-LU-shift 全場相對 L2 分別為 4.62997e-14／4.61774e-14，
KSP iterations 為 124／114；flow 三候選全場相對 L2 皆約 3.33575e-11，
KSP iterations 依 LU／block-ILU／block-LU-shift 為 3／138／121。
這六個作業最大 peak rank RSS=99,385,344 bytes；原序列參考與數值／守恆 gate 保留。
Pressure 三候選亦已通過，總計 9 作業、18 rank reports，見 `regular-acceptance.json`
及 `audit_regular.py`。Pressure LU／block-ILU／block-LU-shift 相對 L2 分別為
9.95969e-14／9.96206e-14／9.96123e-14，KSP iterations 為 2／143／123；
assembly 最大 rank 秒為 27.3637／26.5182／26.0157，linear 為
0.240324／0.0527229／0.0692620。九個作業最大 peak rank RSS=100,704,256 bytes。
後續小切割幾何的結果另見 [小切割評估](HPC_04B_SMALL_CUT_PROGRESS.md)。還須彙整全部成功／失敗原因與場 gate，再處理更小切割
元素、其他網格／rank 與排程資源驗收。128-element 方管四候選已通過，見
[加密方管結果](HPC_04B_REFINED_DUCT_PROGRESS.md)；完整 moving regression 在 rigid wall trace gate 失敗，尚待修復；HPC-04 整體保持未完成。
