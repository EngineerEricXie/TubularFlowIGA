# HPC-04B C2 方管候選預條件器

- 狀態：16-element、2-rank 小案例候選驗證通過；HPC-04B／C 已於 2026-09-11
  由 [完成報告](HPC_04BC_COMPLETION_REPORT.md)簽核。
- 日期：2026-09-09（美東）。Native frozen binary 為 `01aa4c8` 的
  `outputs/hpc04/immersed/native-v4-binary`；fixture／harness 在 `f364a5f` 後新增。
  沒有更動 production solver 或資料格式。本階段不以當前 Makefile 的 source identity
  冒充 frozen native binary 的 identity。

新增 `hpc_duct_solver_fixture`，重用既有 `CouplingFixture.hpp` 的 C2 方管與
Phase 9 的五 domain、0D source、1D 分岔及兩個不同 RCR outlet。可指定橫向／軸向
元素與 packed rank 數，固定兩步、dt=0.01；至少一個元素/rank。以原 `Validate`
檢查完成標記、完整 port／edge／step 身分、strong coupling 收斂、0D constitutive
與儲存量、每步全系統體積平衡，不複製或放寬驗收公式。
Validator 的 `OUTPUT` 父目錄必須含對應 graph `simulation_config.json`。

目前案例為橫向 2 × 2、軸向 4，共 16 elements、175 nodes、700 flow DOFs、2 ranks。
No-slip wall、指定 inlet velocity profile 及 outlet pressure traction 保留；出口 traction
提供絕對壓力參考。此案例不涵蓋純速度邊界的 pressure nullspace、immersed 小切割
元素、ghost stabilization 或所有 port constraint 類型。

## 候選與驗收

所有候選均為 FGMRES、rtol=1e-12，其餘 production 數值預設保留：

| 候選 | PC 配置 |
|---|---|
| lu | LU／MUMPS |
| bjacobi | block-Jacobi，各子群 preonly／ILU |
| schur-lu | full Schur，速度與壓力子區塊各 LU／MUMPS，A11 Schur preconditioner |
| schur-gamg | full Schur，速度子區塊 GAMG，壓力 Jacobi，A11 Schur preconditioner |

Fieldsplit 明確指定 block size=4、velocity fields=0,1,2、pressure field=3。
只對速度子區塊使用 GAMG。Pressure Jacobi 為待評估候選，不預設它較好。

原 validator 的 pressure residual ≤1e-6、flow residual ≤1e-10、3D mass imbalance
≤1e-8、per-step global volume balance ≤3e-15 均保留。每個成功候選另比較所有
accepted checkpoint flow／1D field shards，沿用逐值 `1e-12 + 1e-6*abs(reference)`
及整體 L2 `1e-12 + 1e-6*reference_l2`。非零退出、缺測量或比較失敗都記為失敗；
參考 LU 失敗會停止整個 suite，其他失敗候選保留並繼續下一個。

`scripts/hpc_duct_solver_candidates.py` 保存 fixture/binary hashes、指令、選項、原始
validator 輸出、每 rank RSS／wall time、assembly／solver setup／linear solve／
communication／output 分階段 profile。不同 candidate 使用相同幾何與物理配置。

```bash
make -C solvers/coupling hpc_duct_solver_fixture
python3 scripts/hpc_duct_solver_candidates.py \
  --binary outputs/hpc04/immersed/native-v4-binary \
  --transverse 2 --axial 4 --ranks 2 \
  --output-dir NEW_OUTPUT
```

## 證據與限制

本機 PETSc 3.15.5 real64/int32、Open MPI 4.1.2、GCC 11.4；OMP／OpenBLAS=1。
Native 執行 profiling 開啟，沒有更改原 CPU build。證據根目錄
`outputs/hpc04/pde-candidates`；正式 suite 為 `duct-final-v1`。
執行時另有完整 moving regression，因此本次時間僅是功能驗證的觀測，
不能用於無干擾 speedup 或 strong／weak scaling 結論。

前置 `probe-v1` 使用原一元素 bifurcation fixture（不是加密方管）：四個候選均退出 0
且場比較通過，block-Jacobi／Schur-GAMG 的 history 比較停止於 roundoff 級的
normalized residual 差異。直接比較 normalized residual 的相對誤差沒有代表原守恆
門檻；正式方管 suite 改為呼叫既有 Phase 9 validator。前置 probe 不作正式通過紀錄。

正式 suite 四個作業、8 份 rank report 全部通過；各候選均完成 42 次線性求解。
`audit.json` 另核對兩個 accepted step 的有效 prefix、KSP／PC 與正收斂原因，
保存 fixture／harness source hashes 與 acceptance hash。

| 候選 | 最大場相對 L2 | 累計 linear iterations | assembly s | setup s | linear solve s | peak rank RSS bytes |
|---|---:|---:|---:|---:|---:|---:|
| lu | 參考 | 42 | 30.1305 | 0.913038 | 0.0266137 | 59,211,776 |
| bjacobi | 1.39369e-11 | 4,274 | 30.8797 | 0.0834245 | 0.431851 | 46,645,248 |
| schur-lu | 5.52741e-12 | 62,351 | 31.3199 | 0.0413609 | 37.2117 | 60,715,008 |
| schur-gamg | 1.08564e-10 | 12,776 | 31.3542 | 0.0421752 | 5.69065 | 51,568,640 |

時間為各 rank 的最大 exclusive phase time，分別取最大值，不能相加冒充同一 rank
wall time。Communication 最大值依表序為 0.531543／0.463317／0.886034／0.726576 s，
output 為 0.127022／0.130280／0.140184／0.149620 s。未涉及 CUDA allocation。
四個候選最大 pressure residual 均約 8.13977e-7，最大 flow residual 2.968e-14，
最大 3D mass imbalance 6.25272e-12。其餘 0D／全系統 gate 由原 validator 檢查通過。

Schur／LU 的迭代成本在此案例高，保留為候選的負面證據；不因場精度通過而稱其
有效率。迭代數不受背景 wall-time 干擾，但本次時間仍不足以作一般效能排名。
後續需評估 A11 Schur 近似的限制、網格加密與 rank sweep，評估失敗及負收益配置，
再擴展到其他物理／幾何條件與具排程資源的大案例。維持所有既有求解策略預設。
