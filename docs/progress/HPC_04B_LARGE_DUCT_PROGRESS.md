# HPC-04B 1024 元素方管候選驗收

狀態：LU／block-Jacobi 兩候選通過；HPC-04B／C 已於 2026-09-11 由
[完成報告](HPC_04BC_COMPLETION_REPORT.md)簽核。
日期：2026-09-10。Frozen native binary 對應 `01aa4c8`，啟動時 repository
為 `c920d9e`。此批沒有使用後續 fitted quadrature 或 PVTU 修改。

沿用五 domain 強耦合 C2 方管、0D／1D／RCR、兩步 dt=0.01 與原 Phase 9
物理門檻；網格 8×8×16，共 1024 元素、2299 控制點、9196 flow DOFs，4 ranks。
FGMRES rtol=1e-12，比較 LU／MUMPS 與 block-Jacobi／sub ILU。
兩作業、八份 rank reports 全部 exit 0、無 timeout；accepted step 的有效
solver prefix、PC 與正 convergence reason 均核對。

原 validator 重新執行通過，包含 pressure residual≤1e-6、flow residual≤1e-10、
3D mass imbalance≤1e-8、per-step global volume balance≤3e-15 及其餘
constitutive／coverage 條件。最大 pressure residual=6.26674e-7、flow
residual=2.99986e-16、3D mass imbalance=7.95338e-12。兩步的 1D 與四份
flow shards 重算比較，block-Jacobi 最大 relative L2=3.68215e-11，未放寬門檻。

| 候選 | 累計線性迭代 | assembly s | setup s | linear solve s | peak rank RSS bytes |
|---|---:|---:|---:|---:|---:|
| LU／MUMPS | 42 | 2128.771 | 171.655 | 0.723528 | 412999680 |
| block-Jacobi／ILU | 3012 | 2139.524 | 1.76630 | 11.9201 | 131657728 |

每個候選有 42 次線性 solve。時間為各 phase 最大 rank exclusive seconds，
不能相加當成某 rank wall time。Communication 分別為 32.3034／36.1928 s；
output 為 0.182190／0.181142 s。Wrapper wall 最大值為 2318.256／2171.397 s。
各 rank 峰值 RSS bytes：

- LU：338083840、412999680、327602176、309329920。
- block-Jacobi：122736640、131497984、131657728、122605568。

硬體為本機 i9-14900KF、可見 16 logical CPUs，PETSc 3.15.5 real64/int32、
Open MPI 4.1.2、GCC 11.4、OMP／OpenBLAS=1。無 CUDA allocation。
同時有 moving 與其他回歸，未隔離 CPU affinity；本批提供數值與 RSS 觀察，
不是無干擾 scaling／speedup 證據。較低 RSS 與 setup 成本值得後續固定配置
量測，但不據此變更預設 solver；亦未在此網格評估 Schur 候選。

## 證據與重現

`outputs/hpc04/pde-candidates/large-two-v1/acceptance.json` 保存完整命令、
inputs／binary hashes、rank reports、solver 配置、物理與場驗證。
`launch-source.json` 是保留不修改的啟動紀錄，其 running 欄位不是目前狀態。
完成狀態以 acceptance 與 `audit.json` 為準。`audit_duct.py` 再次檢查所有
input、binary、harness／dependency、stdout／stderr hashes，重跑原 validator
及 field comparator，核對 profile 與每個 rank 的報告。

```bash
python3 scripts/hpc_duct_solver_candidates.py \
  --binary outputs/hpc04/immersed/native-v4-binary \
  --transverse 8 --axial 16 --ranks 4 --timeout 7200 \
  --candidates bjacobi --output-dir NEW_OUTPUT
```

後續完成報告已彙整 fixed-mesh rank sweep、pressure reference／gauge、immersed
條件與候選負收益。跨節點配置仍由 HPC-09 驗收。
