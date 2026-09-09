# HPC-05C：完整 native graph 的新作業續跑

日期：2026-09-09。基準 revision：`9de6ea7a0dd7ca2607cb4085f1bcfb2b67f371d3`
加本批修改。**native graph 的單機整合已通過；05C 保留部分完成，尚待大型與跨節點
排程驗收。** 整份 TODO 的 38 項範圍維持不變。

0D／1D／貼體 3D providers、完整 accepted history、CLI 與全群恢復已接入正式
`iga_multidomain_flow`／`iga_1d_3d_bifurcation`。不中斷與新 MPI 作業續跑的所有
輸出檔案及 checkpoint payload 精確一致，包含每一步的迭代、場、時鐘、pressure
guess、species donor 與帳目。操作說明見 [graph restart](../COUPLED_RESTART.md)。

## 實作

[NativeGraphCheckpoint.hpp](../../solvers/cpu/include/NativeGraphCheckpoint.hpp)
建立完整 domain catalog；0D／1D／coupling controls／history 歸 group root，3D owned
fields 由各 rank 寫自己的分片。receipt 只交換 metadata，不 gather 大場。保存點在
所有 finalize、donor、accepted result、下一步 pressure guess 與時鐘發布之後。

每個 graph 的身分包含全部 configuration text、外部資產 SHA、數值 gates、有效
PETSc options、compiler／OpenMP／fast-math、MPI／PETSc 版本、scalar／index 寬度、
thread 環境及 group local→world rank 映射。絕對 case 路徑不納入內容身分。
[來源指紋工具](../../scripts/native_checkpoint_source_identity.py) 在編譯時 hash native
runner、相關共用 header 集合、Makefile 與自身內容，包含未 commit 修改且不需 `.git`。
不同來源會拒絕載入；docs／測試產物改動不改變此指紋。final source digest 為
`5a2958ffdad9d894fd460eaa5694880887abdffb6e04d6b6760ccedc41cf041f`。

loader 先驗證完整世代，建立未發布的 fresh graph 候選，各 rank 校驗 shared／owned
分片；history 的 domain／port／edge／species catalog、完整 step clock 與最後 pressure／
donor controls 也須相符。local candidate read 全群同意後才執行 typed restore；所有
owner 恢復完成前不能 advance 或發布輸出。恢復中 PETSc 失敗時，整份候選 graph
丟棄，不聲稱可對運行中的 graph 做原子替換。

新建 checkpoint root 的各層目錄與既有 parent 都同步。epoch 遇到既有完整或不完整
路徑時選新 ID，最多 65,536 次，絕不覆蓋原世代。每個分片先完成 checksum 與 sync，
root 驗證完整 catalog 後最後發布 manifest。故意寫入錯誤的 hooks 僅存在於
`IGA_COUPLED_CHECKPOINT_TESTING` 測試建置，正式 runner 沒有新的故障注入開關。

[GraphAcceptedHistory.hpp](../../include/GraphAcceptedHistory.hpp) 的兩種 version-1
records 保存 flow／species result 的所有成員，包含完整 hydraulic iterations、port
optionals、原 graph 方向的 donors、canonical diagnostics、transport order、domain／global
species accounting 與 0D state／accounting。species diagnostics 排序 port refs 時會反向
標示 donor；loader 保留這個與 runtime donor map 不同的方向語義。
[CheckpointRecordStream.hpp](../../include/CheckpointRecordStream.hpp) 以 uint64 little-endian
長度 framing 逐筆寫入，每筆 metadata 上限 16 MiB；完整 history 可超過此值。

[訊號處理](../../solvers/cpu/include/NativeGraphCheckpointSignal.hpp) 從輸入前安裝
SIGUSR1 handler，只設定 `sig_atomic_t`。任一 rank 的 request 在下一個 accepted
boundary 全群傳播，強制保存後正常停止；預警不需要碰巧落在 `--checkpoint-every`。
正常最後步與 `--stop-after-step` 指定的最後步也保存。退出後恢復先前 signal action。

## 驗收門檻與結果

門檻在比較前固定：相同分區的完整輸出 bytes、每代所有 shard／manifest bytes 精確
一致；中斷重試產生新 epoch ID 時，逐片 payload 精確一致。時間、IDs、presence、
counts 及 history prefix 必須完整，numeric CSV 必須有限。原有 Newton、mass、species
與 coupling convergence gates 不變。場相對 L2 為 0，因全部 FP64 field bytes 一致。
不同 rank 數各自與其不中斷 reference 比較，不要求跨 rank reduction 逐位元一致。

| 正式 v6 驗收 | 結果 |
|---|---|
| 三 ranks 完整 suite | 63 個新 MPI 作業：41 成功、22 預期失敗；成功作業 123 份 rank reports，exit 0、無 timeout、stderr 空 |
| 七種圖／方法 | 0D–3D、1D–3D flow、兩物種、0D fixed relaxation、1D Aitken、近零 donor 後反轉、兩個 3D islands；full／save-3／resume 的完整輸出與 bundle bytes 精確一致 |
| 方法辨識 | fixed／Aitken 有多次 iteration；第 4 步第一個 applied pressure 不等於初始 guess，恢復後歷史仍完全一致 |
| donor 辨識 | 第 3 步非零 donor 在第 4 步近零流量保留，第 5 步反轉；permuted graph 的 canonical diagnostics 正確 |
| 三個中斷點 | zero／flow／species 各測 precommit exit 1、分片寫入 5 bytes 後 exit 86、manifest.tmp 完成但發布前 exit 87；各以新作業恢復，舊完整世代 hashes 不變 |
| rank-local 寫入失敗 | rank 1 寫自己的 flow 分片前拋例外，全群 exit 1；新的作業仍恢復上一完整世代 |
| 最新世代不可用 | byte corruption、缺片、跨 epoch 換片；由第 2 步恢復，再生成第 3 步的新 ID，原錯誤世代內容不變 |
| 不相容拒絕 | source SHA、Newton controls、PETSc options、case、rank 數、field 順序、外部 waveform table，以及短於 checkpoint 的 horizon；無新 graph 輸出 |
| 使用錯誤 | interval 沒有保存目錄、空 checkpoint root 都 exit 1 |
| 搬移／唯讀 | 移動 case 仍恢復；不指定保存目錄時原 bundle hashes 不變；僅輸出 manifest 的絕對 case 路徑允許不同 |
| 排程預警機制 | 只對 root 發 SIGUSR1：第 3 步保存中收到則第 4 步保存停止；輸入載入中收到則第 1 步保存停止；頻率皆為 3，新作業續跑精確一致 |
| 1／2／4 ranks | 各跑 zero／flow／species／multi，12 作業，共 36 作業、84 份 reports；每個配置均通過 full／save／resume |
| 獨立 communicator | 三個新 WORLD-3 作業，分成 1＋2 ranks，分別跑 zero 與 species 並使用不同 PETSc option database；兩群各自完整輸出與 bundle 精確一致 |
| 既有五域循環 | Phase 9 source 0D→3D→分叉 1D→兩個 RCR 0D，在 1／2 ranks 各跑 full／save-3／resume，沿用原 fixed controls，八步全部輸出與分片 bytes 精確一致 |
| history codec | 3,128 個預期拒絕；21,168,000-byte 完整 history 超過 16 MiB 仍精確 roundtrip，該測試最大 chunk 1,315 bytes |
| 儲存回歸 | 原 bundle 測試 734 個拒絕與 32 MiB 串流仍通過 |

`wrong-source`／`wrong-field-order`／`wrong-external-table` 明確於 checkpoint discovery
拒絕；wrong-ranks 在 `.ntiga` domain preflight 拒絕，沒有錯誤分區恢復。

final-source native ASan／UBSan 的 zero／flow／species full／save／resume 共 9 作業、
27 份 rank reports 全通過且 stderr 空；所有輸出與分片 bytes 精確一致。bifurcation
相容入口另跑三個新作業也通過。checkpoint 關閉時，final native binary 的四個模式
與封存舊 executable 的既有相同 case 輸出逐 byte 一致（`legacy-v6` 對 `legacy`）。
單元 history ASan／UBSan／LSan 亦通過；native MPI sanitizer 的 leak detection 關閉，
因系統 MPI／PETSc 未以 sanitizer 編譯，不宣稱 third-party leak audit。

## 重現與環境

本機 Linux x86_64、GCC 11.4／OpenMPI 4.1.2／PETSc 3.15.5 real64/int32，
`-O3 -std=c++17 -Wall -Wextra -Wpedantic -fopenmp`；每 rank 一 OpenMP／BLAS thread，
preonly／LU／MUMPS，未設定 affinity。C++ 建置無新 warning。受限建置中的
`opal_ifinit errno=1` 為 socket 權限訊息；MPI／sanitizer／strace runtime 使用本機所需權限。

```bash
make -C solvers/coupling iga_multidomain_flow iga_1d_3d_bifurcation \
  native_graph_checkpoint_fault_test native_graph_checkpoint_fixture_test \
  native_graph_checkpoint_groups_test graph_checkpoint_history_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
solvers/coupling/graph_checkpoint_history_test
python3 scripts/hpc_native_graph_checkpoint.py --output-root NEW_SUITE_ROOT
python3 scripts/hpc_native_graph_checkpoint.py --output-root NEW_RANK_ROOT \
  --positive-only --ranks 4 --modes zero,flow,species,multi
```

`NEW_*_ROOT` 必須尚不存在；`--ranks` 可改 1／2／4（完整 fault suite 至少兩 ranks）。
script 由 repository fixture builder 產生 configuration／partition，保存二進位 SHA、
所有 command argv、expected／actual exit code 與 checks。正例經有 timeout 的
`hpc_rank_run.py` 收集各 rank wall／RSS／stdout／stderr hashes。所有預期非零作業
直接使用有 timeout 的 mpiexec，避免 OpenMPI 終止 wrapper 後來不及寫 `run.json`。

證據根目錄為 `outputs/hpc05/native-graph/`：`final-suite-v6`、`ranks-{1,2,4}-v6`、
`groups-v6`、`io-v6`、`sanitized-suite-v6`、`bifurcation-v6`、`circulation-{1,2}`、
`legacy-v6`，以及 history／bundle／build logs。五域案例由封存的
`circulation-fixture.cpp` 呼叫既有 Phase 9 `PrepareCase(ROOT, 0.05, 8, RANKS, "root.ntiga")`
生成；原 test 的模型、兩個 C2 duct elements 與 fixed controls 沿用。`run-groups-v6.py` 準備
`group-0/fixture`（zero，1 rank）與 `group-1/fixture`（species，2 ranks），再用同一
`native_graph_checkpoint_groups_test ROOT full|save|resume` 跑三個獨立 WORLD-3 jobs。
每群使用不同 `-ksp_rtol`，確保不意外比較 WORLD options。所有生成資料與二進位
只存 ignored outputs，未提交。

最終驗收合計 126 個作業：104 成功、22 個預期失敗；成功作業共 279 份 rank reports。
`acceptance.json` 彙整各 suite 結果、215 份來源 SHA 與 9,464 份證據 SHA；另有
commit 後的 revision／clean-worktree 核對。早期預覽與失敗 suite 保留原狀，沒有改寫
為成功紀錄。

## 量測與限制

`IGA_PROFILE=1` 的 rank-local profiles 分開記錄 exclusive assembly、solver setup、
linear solve、明確埋點的 application communication 及 output。checkpoint Save／Restore
納入 output；MPI／PETSc library 內部通訊仍留在呼叫者 phase，不是完整網路通訊量。

`run-io-v6.py` 的兩個三 rank 作業，以 `strace -ff -yy -T -e trace=read,write,fsync`
按 rank 記錄 checkpoint 路徑的 syscall，另計新 root 的 parent fsync；traced 續跑的
完整輸出及 bundle bytes 仍與正式 reference 精確一致。`measurements.json` 記錄每 rank
profile／wall／peak RSS 與 checkpoint read／write／fsync 的 calls、bytes、seconds。
這是 syscall latency，未包含 encoding、SHA、其他 metadata syscall；trace 本身有成本，
不能當作未追蹤的 I/O 吞吐，更不能用來宣稱加速。

| traced 作業／rank | assembly s | setup s | linear solve s | communication s | read／write syscall s | fsync s | peak RSS bytes |
|---|---:|---:|---:|---:|---:|---:|---:|
| resume/rank-0 | 0.178781 | 0.018769 | 0.001228 | 0.024182 | 0.535859／0.005547 | 0.308831 | 44,806,144 |
| resume/rank-1 | 0.169060 | 0.018777 | 0.001228 | 1.875285 | 0.077375／0.000421 | 0.033438 | 45,731,840 |
| resume/rank-2 | 0.176998 | 0.018772 | 0.001228 | 1.868434 | 0.077536／0.000629 | 0.032557 | 47,292,416 |
| save/rank-0 | 0.138783 | 0.011029 | 0.001078 | 0.014834 | 0.343704／0.004296 | 0.333398 | 44,937,216 |
| save/rank-1 | 0.141457 | 0.011029 | 0.001075 | 1.520104 | 0.000000／0.000465 | 0.036436 | 45,953,024 |
| save/rank-2 | 0.137797 | 0.011036 | 0.001078 | 1.524701 | 0.000000／0.000416 | 0.037783 | 46,850,048 |

第一代三 rank flow／species fixture 每個 3D island 只有一 cubic element／64 nodes。
1／2／4-rank 和 split 是 ownership／一致失敗／恢復正確性測試，不是 scaling 證據。
部分作業與建置重疊，端到端 wall 也不作效能比較。CUDA allocation 為 N/A。
既有 replicated history vectors、1D state 與 3D boundary arrays 仍存在，恢復另需要
未發布候選；每世代保存完整 prefix 的空間成長與大量小檔成本由 HPC-06 接續。

開發期間曾把 species canonical diagnostics 當成 graph donor 方向而拒絕正常恢復，
已修正並由 permuted／反轉案例驗證。另一個 1D flow fixed-relaxation fixture 首步
MAX_NEWTON 失敗，封存舊 native binary 也重現相同失敗；沒有放寬收斂門檻，fixed
驗收使用已有效的 0D graph。第一版 suite 曾因 expected-exit 作業的 wrapper report
未生成而失敗，故負例改用直接 launcher；失敗紀錄保留在 preview／final-suite。

大型 graph、實際跨節點 checkpoint filesystem 與 scheduler signal forwarding 尚未
執行。依狀態契約的「規模與資源」gate，05C 保持未勾選；下一步使用配置好的 compute
allocation 做這些驗收，不能把本地 MPI 當作跨節點證據。immersed／moving／FSI、VCA
及不同 rank 數恢復仍由 05D 接續，HPC-01C／01D 與其他未完成工作亦維持原範圍。
