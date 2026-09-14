# HPC-05C／07C／07D／09A–D：Bridges-2 跨節點驗收

- 狀態：passed；完成日期（UTC）：2026-09-14。
- 本批七項剩餘驗收完成；原有工作站與 GPU 紀錄沿用各自歷史報告，本批沒有重跑 GPU。
- 乾淨數值驗收 source：`f549be47eb4552119508627b3b2eca0d0ebe9660`；起始乾淨 source：`faf11a3ad2a0dd3197ab662a4ab854a58218f335`。
- Branch：`nextgen-phase-0-runtime-foundation`；account：mch260002p。
- Evidence root：`/ocean/projects/mch260002p/thsieh1/TubularFlowIGA-hpc-acceptance`；最終判定：`cross-node-acceptance.json`。

## 環境與排程

- RM nodes，AMD EPYC 7742，128 physical cores/node、兩個 NUMA domains；實測 hardware 在 graph-v2/attempt-1/hardware-*.txt。
- Open MPI 4.0.5／GCC 10.2.0；PETSc 3.22.1 real64／32-bit indices，MUMPS 5.7.3、ScaLAPACK 2.2。
- HDF5 1.12.1 serial C library，來自 anaconda3-2024.10-1；MPI wrappers 與 PETSc 使用相同 GCC/Open MPI ABI。
- PETSC_DIR=`/ocean/projects/mch260002p/thsieh1/TubularFlowIGA-hpc-acceptance/petsc`；PETSC_ARCH=`arch-bridges2-acceptance`。
- Compute allocation 45782671 內完成 PETSc 與 make hpc-cross-node-binaries；沒有在 login node 執行 simulation。
- ABI／linkage／binary hashes：build/binaries-v2.json；建置輸出：build/partition-fix/final-build.stdout／stderr。
- 每 rank 一 core，OMP／BLAS threads=1；Graph 與最大 scaling 使用每節點 128 ranks。
- Graph → FSI → scaling → finalizer 為 afterok 序列，大型作業最多同時占用兩節點。

| 階段 | Job ID | 實際配置 |
|---|---|---|
| Graph | 45785484 | 7 attempts；256 ranks、每次兩節點 |
| FSI | 45831818 | writer: r285.ib.bridges2.psc.edu, r414.ib.bridges2.psc.edu；reader: r285.ib.bridges2.psc.edu, r414.ib.bridges2.psc.edu |
| Scaling | 45877533 | r401.ib.bridges2.psc.edu, r440.ib.bridges2.psc.edu；ranks 1/64/128/256 |
| Finalizer | 45877534 | 1 CPU，三個正式作業成功後執行 |

## Graph：訊號、checkpoint、requeue 與 restart

- 16,384 elements、24,187 nodes、5 domains、8 macro-steps；dt=0.01。
- 初次受控 SIGUSR1 在第一步接受後發布 checkpoint，Slurm requeue 後在另一組兩節點恢復。
- 各 attempt 的 nodes.txt、build.json、scheduler.json、solver.stdout／stderr 與完整共享 checkpoint 均保留。
- MPI binding 逐 rank 核對七次 attempts：每次 256 ranks、每節點各 128；詳見 graph-v2/mpi-binding-audit.json。
- 初始 staging 逐檔雜湊一致；step 1 checkpoint 的原始 bytes 在後續 publication 後仍不變。
- 單一步驟可超過五分鐘；後續提前送 SIGUSR1，等待 accepted-step checkpoint，沒有在 trial 中強制發布。
- 此次 source 的 launcher 重設 requeue flag；實跑於 launcher 啟動後重新啟用。最終腳本修正將相同操作放在 scontrol requeue 前。
- 最大 normalized pressure residual=`9.28324e-07`（門檻 1e-6）；flow residual=`9.96736e-11`（門檻 1e-10）。

| Step | Coupling trials | 3D mass imbalance（m³/s） | External outward flow（m³/s） |
|---:|---:|---:|---:|
| 1 | 21 | 1.96975e-13 | -6.22853e-14 |
| 2 | 31 | 1.05076e-12 | -1.48991e-13 |
| 3 | 26 | 3.3586e-13 | -2.1811e-13 |
| 4 | 20 | -4.93093e-13 | -3.48931e-14 |
| 5 | 23 | 7.5395e-13 | 2.83682e-14 |
| 6 | 21 | -4.03436e-13 | 9.37427e-14 |
| 7 | 20 | 3.11378e-13 | 2.33375e-14 |
| 8 | 21 | 5.96372e-13 | -4.11314e-14 |

| Attempt | Hosts | Scheduler status |
|---:|---|---|
| attempt-0 | r297.ib.bridges2.psc.edu, r373.ib.bridges2.psc.edu | checkpointed |
| attempt-1 | r195.ib.bridges2.psc.edu, r251.ib.bridges2.psc.edu | checkpointed |
| attempt-2 | r087.ib.bridges2.psc.edu, r472.ib.bridges2.psc.edu | checkpointed |
| attempt-3 | r121.ib.bridges2.psc.edu, r407.ib.bridges2.psc.edu | checkpointed |
| attempt-4 | r409.ib.bridges2.psc.edu, r473.ib.bridges2.psc.edu | checkpointed |
| attempt-5 | r319.ib.bridges2.psc.edu, r419.ib.bridges2.psc.edu | checkpointed |
| attempt-6 | r319.ib.bridges2.psc.edu, r339.ib.bridges2.psc.edu | passed |

## FSI：two-node scheduled tier 與 4→2 paired restart

- Scheduled tier 實際 passed，沒有 skipped；nonzero strong writer 四 ranks 跨兩節點，reader 兩 ranks 跨兩節點。
- Checkers 驗證 fields、surface、traction、Aitken／完整 strong histories、ports、conservation 與 input checkpoint bundle 不變。
- 結構矩陣仍採 bounded single-owner solve；沒有宣稱全分散膜矩陣或跨節點加速。

| 比較 | 最大 field scaled L2 | 最大 history error / tolerance |
|---|---:|---:|
| reference-1 | 0 | 0 |
| writer-4 | 1.23017e-12 | 0.000219258 |

4→2 reader：accepted scaled L2=`7.21654e-13`；next-step scaled L2=`8.95107e-13`；continuation iterations=7。

## Formal strong／weak scaling

- Strong 固定 16,384 elements；weak 固定 256 elements/rank；四種 rank 數各至少三次 fresh process runs。
- 每次 geometry／native physical validation 通過；strong fields 對一 rank reference 通過原訂 relative L2 1e-6／zero-reference absolute L2 1e-12。
- Wall 為每次最大 rank process wall 的中位數與範圍；RSS 為每次最大 rank peak RSS 的中位數。全部原始樣本與各 phase 在 scaling-v6/summary.json。
- Strong efficiency=T1/(p×Tp)；weak efficiency=T1/Tp。Communication 是最大 rank exclusive phase 的中位數，不能把各 phase 最大值相加當作 wall time。

| Mode | Ranks | Elements | Wall median [min,max] s | Speedup | Efficiency | RSS MiB | Iterations | Communication s |
|---|---:|---:|---|---:|---:|---:|---|---:|
| strong | 1 | 16384 | 4689.21 [4678,4689.83] | 1 | 100% | 3902.17 | 147/147/147 | 0.887934 |
| strong | 64 | 16384 | 729.946 [724.992,730.688] | 6.42405 | 10.0376% | 305.68 | 465/465/465 | 602.602 |
| strong | 128 | 16384 | 494.997 [494.949,499.051] | 9.47322 | 7.40095% | 293.043 | 592/592/592 | 427.537 |
| strong | 256 | 16384 | 381.903 [381.491,392.309] | 12.2785 | 4.7963% | 299.055 | 784/784/784 | 346.18 |
| weak | 1 | 256 | 73.7993 [73.5982,74.3988] | 1 | 100% | 215.66 | 31/31/31 | 0.00247269 |
| weak | 64 | 16384 | 725.409 [724.662,729.649] | 0.101735 | 10.1735% | 308.117 | 465/465/465 | 601.739 |
| weak | 128 | 32768 | 995.528 [994.233,997.157] | 0.0741308 | 7.41308% | 503.805 | 1385/1385/1385 | 856.691 |
| weak | 256 | 65536 | 688.724 [686.214,702.145] | 0.107154 | 10.7154% | 900.383 | 549/549/549 | 574.364 |

Strong velocity 最大非零 reference relative L2：`3.75761e-14`。

Strong pressure 最大非零 reference relative L2：`6.64998e-13`。

## 修正與失敗證據

- Graph job 45784010 的 fixture index 在每 rank 列出所有元素，診斷後取消。改為 production packer 的 owned-row required-element 規則；physical records 不變，索引 oracle 與四 rank 場 byte comparison 通過。
- Open MPI 4.0.5 改用 ess_base_forward_signals SIGUSR1；修正 spool repo 路徑、PSC LOCAL staging、validator build target 與 weak efficiency。
- 初始 PETSc 缺少 MUMPS；於 compute allocation 另建相容安裝。Requeue ordering 另有完整 wrapper regression。
- 診斷保留於 graph-v1/、preflight/、build/partition-fix/、final-fix/ 與 slurm/；沒有放寬數值門檻或以工具小案例替代 formal scaling。
- Installed Slurm 的 sstat 不支援 Elapsed／CPUTime；原始 probe 錯誤保留，另存 sacct JSON；最終 recorder 改用已實測的 AveCPU 欄位。
- Finalizer 首次 2G/core 申請超過 PSC 2000 MB/core，submission 失敗且未取得 job ID；1900 MB 重提成功，兩份 submission JSON 保留。
- FSI 首次 job 45793958 在三小時 allocation 上限於 reader 階段 TIMEOUT；reference 與 writer 已通過 strong checker。保留 fsi-v1/、slurm/fsi-45793958.* 及 terminal-evidence/45813868/，以相同乾淨 source 完整重跑，僅將 allocation 改為 06:30:00；每個 process 7200 秒與數值門檻不變。
- 本報告如實保留低效率配置；binary、database、case、results、scheduler logs 與 JSON evidence 留在 shared root。
- Scaling 啟動失敗 45831824／45835451：Bash case 預設值未加引號，產生損壞的 rank pattern；環境覆寫仍留下未配對括號。兩次均在正式 solver 前失敗，stdout/stderr、scheduler JSON 與診斷保存在 scaling-v2.*、scaling-v3.* 及 slurm/。提交 wrapper 修正兩個預設值引號，並加入預設／自訂路徑回歸測試。
- Scaling 45835866 在六個樣本通過後，完整 weak128（32,768 elements）線性求解於 5,000 次迭代仍未收斂；scaling-v4/failure-diagnosis.json 與原始紀錄保留。ASM overlap2／ILU(1) 診斷 45839794 觸發 Slurm OOM；此前 45839653 的 viewer YAML 旗標解析錯誤亦已保存與修正。
- 完整 weak128／weak256 診斷 45839897 通過：保留原 block-Jacobi／ILU(0)，僅將 FGMRES restart 設為 200，原線性 rtol 1e-8、5,000 次上限、Newton 與守恆門檻不變。正式三輪完整重跑採同一設定，透過 PETSC_OPTIONS_YAML 配置；實際 KSP view 在每個 rank stdout。診斷證據位於 scaling-diagnostic-gmres200-45839897/。
- Scaling 45841321 因 Slurm NODE_FAIL 在首個樣本中止；沒有完成樣本可計入驗收，rank stderr 沒有數值錯誤。scaling-v5/failure-accounting.json、failure-diagnosis.json 與原始 stdout/stderr 完整保留；重新分配節點後，以相同設定完整重跑。
- 正式 scaling 提交副本及雜湊封存於 scaling-v6/submission-wrapper.sbatch 與 wrapper-provenance.json；與 f549 source 的差異為 Bash 預設值引號及上述 restart 環境配置。solver、checkers、build 維持乾淨 f549 revision，完整 cases/scaling-v2 的物理輸入不變。這些 wrapper 修正於完成 commit 納入 source。

重跑 finalizer：

```bash
python3 scripts/hpc_finalize_cross_node.py \
  --graph-output /ocean/projects/mch260002p/thsieh1/TubularFlowIGA-hpc-acceptance/graph-v2 \
  --fsi-output /ocean/projects/mch260002p/thsieh1/TubularFlowIGA-hpc-acceptance/fsi-v2 \
  --scaling-output /ocean/projects/mch260002p/thsieh1/TubularFlowIGA-hpc-acceptance/scaling-v6 \
  --output /ocean/projects/mch260002p/thsieh1/TubularFlowIGA-hpc-acceptance/cross-node-acceptance.json
```
