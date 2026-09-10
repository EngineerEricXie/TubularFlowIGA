# HPC-04B 加密方管候選結果

- 狀態：128-element、4-rank 的四候選驗收通過；HPC-04 整體仍進行中。
- 日期：2026-09-09（美東）。Native frozen binary 為 `01aa4c8`，fixture builder 為
  `42d4d39` 的版本；此批執行的是之後擴充候選、但尚未加入 stderr 保存的 interim
  harness。`harness-at-run.py` 已依記錄的 SHA256 精確還原保存，不能歸為目前
  harness revision 的執行結果。

使用原五 domain 強耦合 C2 方管，橫向 4×4、軸向 8，共 128 elements、539 nodes、
2156 flow DOFs、4 ranks，兩步 dt=0.01。原 pressure traction reference、0D／1D／
RCR、非線性與物理 gates 保留。每個候選與同網格 LU 參考比較 checkpoint field。
此處網格與 rank 都相較先前 16-element／2-rank 作業改變，因此不是隔離單一變因的
mesh convergence 或 scaling 研究。

四個作業、16 份 rank report 均退出 0；每個作業有 42 次線性 solve。
原 Phase 9 validator 的 pressure residual≤1e-6、flow residual≤1e-10、3D mass
imbalance≤1e-8、per-step global volume balance≤3e-15 及其餘 constitutive／coverage
條件全部通過。最大 pressure residual=6.57641e-7、flow residual=1.00479e-14、
3D mass imbalance=1.78925e-11。

| 候選 | 最大場相對 L2 | 累計線性迭代 | assembly s | setup s | linear solve s | peak rank RSS bytes |
|---|---:|---:|---:|---:|---:|---:|
| LU／MUMPS | 參考 | 42 | 262.400 | 11.5438 | 0.107621 | 106,205,184 |
| block-Jacobi／ILU | 2.97464e-11 | 3,724 | 261.341 | 0.120097 | 1.37144 | 54,841,344 |
| A11 Schur／LU | 7.37924e-11 | 84,362 | 260.141 | 0.149266 | 264.737 | 93,114,368 |
| A11 Schur／velocity GAMG、pressure Jacobi | 1.67480e-10 | 10,553 | 282.851 | 0.166418 | 20.6464 | 75,128,832 |

時間為各 phase 的最大 rank exclusive seconds，不能相加冒充同一 rank wall time。
Communication 最大值依序 6.91567／8.79434／8.16911／7.51697 s；output 為
0.164843／0.142199／0.164383／0.256278 s。工作站可見 16 logical CPUs（Intel
Core i9-14900KF）、PETSc 3.15.5 real64/int32、Open MPI 4.1.2、GCC 11.4，
OMP／OpenBLAS=1。執行時有其他回歸作業，僅作功能驗收觀測，不宣稱 speedup。
未涉及 CUDA allocation。較高 Schur／LU 迭代成本仍保留為負面結果，未改 production
預設或放寬停止／場門檻。

## 證據與重現

本機 `outputs/hpc04/pde-candidates/refined-four-v1/acceptance.json` 保存各 job 指令、
binary／fixture hashes、原 validator 輸出、場比較及 profile；`audit.json` 檢查
accepted step 有效 KSP／PC／prefix／正 reason，並保存 execution harness hash 與
acceptance hash。`audit_duct.py` 可重新核對彙整。

```bash
python3 scripts/hpc_duct_solver_candidates.py \
  --binary outputs/hpc04/immersed/native-v4-binary \
  --transverse 4 --axial 8 --ranks 4 --timeout 1200 \
  --output-dir NEW_OUTPUT
```

目前 harness 保留同一四候選預設，並額外保存失敗 diagnostics。上述命令重跑將是
新 harness 的新證據，不覆寫本批結果。

還須固定網格的無干擾 rank sweep、其他 pressure reference／nullspace、immersed
切割條件、更大網格與排程資源驗收。不能由這兩種網格的候選結果宣稱 HPC-04B／C 完成。

後續 1024-element、4-rank 的 LU／block-Jacobi 已完成並通過原門檻，見
[較大方管進度](HPC_04B_LARGE_DUCT_PROGRESS.md)。
