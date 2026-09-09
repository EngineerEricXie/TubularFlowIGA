# HPC-05C：貼體 3D 分片與新 MPI 作業恢復

日期：2026-09-09。基準 revision：`91c06c4117922c29d1cc6e75acf16dfe1d651693`
加本批修改。狀態：**部分完成，HPC-05C 保持未勾選**。

貼體 flow／transport 已能把 accepted state 寫入 bundle 分片，終止原作業後，
由新的三 rank MPI 作業恢復第 3 步並精確接續第 4–6 步。穩態、Backward Euler、
Backward Euler＋RC outlet 三種模式通過；寫場分片中或 manifest 發布前終止，
都不會改變上一份完整 checkpoint。此批完成 runtime 儲存元件及 MPI 驗收，
native graph 的 provider／歷史／CLI 尚未接入。

## 實作與格式

[CheckpointWordStream.hpp](../../include/CheckpointWordStream.hpp) 提供 little-endian
uint64／finite IEEE binary64 串流。writer buffer 為 64 KiB；reader 可接收任意 byte
fragment，只保留未滿一個 word 的 bytes，不累積整份 payload。計數運算檢查溢位，
拒絕 trailing／truncated bytes 與非有限場值；Consume 或 sink 失敗後禁止重用串流。

[BodyFittedAcceptedCheckpoint.hpp](../../solvers/cpu/include/BodyFittedAcceptedCheckpoint.hpp)
將 metadata 與大陣列分開：

| 分片 format | 內容與範圍 |
|---|---|
| `body-flow-v1` | `IGA_BODY_FITTED_FLOW/1` metadata：configuration SHA、accepted clock／macro dt、累積 linear iterations、global rows、boundary counts、tractions、outlet models 與動態值 |
| `body-transport-v1` | `IGA_BODY_FITTED_TRANSPORT/1` metadata：configuration SHA、真實 accepted steps、dt、global rows |
| `body-boundaries-v1` | 兩個 uint64 counts，接 velocity triples 與 pressure values；每 domain 只寫一次 |
| `owned-real-field-v1` | global rows／owned begin／owned end 三個 uint64，接該 rank 的 owned FP64 values；flow、transport 各自一片 |

metadata 延用 16 MiB 上限及 bounded maps。大場不受 metadata 大小上限限制，
仍受 bundle 的 payload 上限約束。immutable boundary masks、scalar boundary 配置
與拓撲從 verified target 建立；typed restore 會再驗證配置與 replicated state。
既有 runtime 內的 replicated boundary arrays 仍存在，HPC-06A 尚須處理其記憶體成本。

兩個 runtime 的 `CreateCheckpointRestoreCandidate()` 只允許 checkpoint-bound、
open、fresh Committed owner，依其可信形狀建立未發布 decoder storage。loader 不依
磁碟提供的 count 配置場大小；先核對 clock／identity／catalog／shape，再逐片填入，
完成 bundle SHA 檢查後才能 typed restore。任何 local I/O 失敗經 group agreement
傳播；已部分填入的 candidate 必須丟棄。這不是 live graph 多個 owner 的原子替換。

[BodyFittedCheckpointBundle.hpp](../../solvers/cpu/include/BodyFittedCheckpointBundle.hpp)
提供 **local** catalog、producer、loader。`world_ranks[i]` 指定 solver local rank i
對應的 bundle rank；shared shards 歸 local rank 0，owned shards 歸對應 world rank。
三 rank flow＋transport 為九片。場不經 MPI gather；
[GatherCheckpointReceipts](../../solvers/cpu/include/CollectiveCheckpointReceipts.hpp)
只彙整 bounded receipt metadata，支援 empty producer，檢查 epoch agreement、總片數
與重複 IDs。root 依完整 catalog 重讀校驗所有分片，再最後發布 manifest。

載入時 root 先 discovery／驗證完整 bundle，廣播 canonical manifest；各 rank 再讀
shared arrays 與自身 owned fields，各分片讀取仍核對 SHA。future graph coordinator
必須提供完整跨 domain catalog、compatibility identity、全作業候選發布與 histories。

## 驗收

[測試程式](../../solvers/cpu/tests/test_body_fitted_checkpoint_bundle.cpp) 重用
accepted-state fixture：一個 cubic element、64 nodes，加 tracer 的 time derivative、
volume source 與 diffusion。比較門檻事先固定為精確一致：metadata／owned fields／
boundaries／outlets／clock／iterations／transport steps 及 `TotalMass()` 都納入 SHA。
每種模式保存第 3 步，保存作業繼續算到第 6 步產生 reference；新的 launcher 從磁碟
恢復，剛恢復及其後三步都逐 rank 比較 reference。這是小型正確性驗收，沒有 scaling 宣稱。

| 檢查 | 結果 |
|---|---|
| 新作業 restart | 三種模式的第 3–6 步完整狀態與 mass 精確一致 |
| rank 映射 | RC 模式 solver local→world map 為 `{2,0,1}`；shared owner 為 world rank 2，恢復與 owned interval 均通過 |
| 單 rank 載入失敗 | 只有 local rank 1 指向 missing root；全群拒絕，fresh runtime counts 保持 0；重建 decoder candidate 後可正常恢復 |
| receipt agreement | empty producers 正常；單 rank epoch 不一致、跨 rank duplicate shard ID 全群拒絕，沒有發布錯誤 manifest |
| 寫分片中終止 | 寫入 flow payload 5 bytes 後 `_exit(86)`；MPI 作業退出 86，留下 `.tmp`，無 manifest |
| manifest 前終止 | 九片完成後 `_exit(87)`；MPI 作業退出 87，九份 `.shard` 存在，無 manifest |
| 中斷後新作業 | 前兩種模式各排除一個不完整世代，回到第 3 步；三種模式仍精確接續到第 6 步 |
| 保護舊世代 | 三份既有 manifest 的 SHA 在兩次中斷前後不變 |
| codec | 910 個預期拒絕，包含 metadata／field／boundary 的所有截斷位置、錯誤 clock／shape／ownership／identity、trailing bytes、Inf 與失敗串流重用 |
| 大／空 owned field | 24-byte 空場與 25,165,848-byte（24 MiB＋header）場 roundtrip 精確一致，最大 chunk 65,536 bytes |
| typed-state 回歸 | 重建原測試；WORLD 3-rank 與獨立 1／2-rank groups、三種模式共九組精確恢復與原有負例通過 |

release 與 ASan／UBSan codec 都通過 910 個拒絕。sanitizer 另以新的三 rank
save／resume 作業驗證相同三種模式；系統 MPI／PETSc 未經 sanitizer 建置，設定
`ASAN_OPTIONS=detect_leaks=0 UBSAN_OPTIONS=halt_on_error=1`，不宣稱 third-party leak audit。

## 命令、環境與證據

Linux x86_64 本機工作站，GCC 11.4／OpenMPI 4.1.2／PETSc 3.15.5 real64/int32，
三 ranks、每 rank 一 thread、preonly／LU／MUMPS，未綁定 affinity。C++17 warnings
建置無新 warning；restricted build 的 `opal_ifinit errno=1` 為 socket 權限訊息，
MPI 與 sanitizer runtime 以所需本機權限執行。

```bash
make -C solvers/cpu body_fitted_checkpoint_bundle_test body_fitted_accepted_checkpoint_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps'
solvers/cpu/body_fitted_checkpoint_bundle_test --codec unused
mpiexec --oversubscribe -np 3 solvers/cpu/body_fitted_checkpoint_bundle_test --save NEW_ROOT
mpiexec --oversubscribe -np 3 solvers/cpu/body_fitted_checkpoint_bundle_test --resume NEW_ROOT
timeout --kill-after=5s 120s mpiexec --oversubscribe -np 3 \
  solvers/cpu/body_fitted_checkpoint_bundle_test --crash-stream NEW_ROOT
# 預期 exit 86；接下來的 crash-manifest 預期 exit 87。
timeout --kill-after=5s 120s mpiexec --oversubscribe -np 3 \
  solvers/cpu/body_fitted_checkpoint_bundle_test --crash-manifest NEW_ROOT
mpiexec --oversubscribe -np 3 solvers/cpu/body_fitted_checkpoint_bundle_test --resume-interrupted NEW_ROOT
```

每組 save 要使用全新的 ROOT；resume 與 interruption 使用同一 ROOT。正式證據的
成功 MPI 作業用 `scripts/hpc_rank_run.py` 收集每 rank stdout／stderr／RSS／wall，
故意退出作業直接由有 timeout 的 mpiexec 啟動，避免 rank wrapper 的獨立 session
影響 launcher 終止語義。sanitizer 使用 `-O1 -g -fsanitize=address,undefined
-fno-omit-frame-pointer`，save／resume 用另一個全新 ROOT。

證據位於 `outputs/hpc05/body-fitted-bundle/`：`save-mpi`、`resume-mpi`、
`recovery-mpi`、`typed-mpi`、`sanitized-save-mpi`、`sanitized-resume-mpi`，另有
codec logs、兩個 crash logs、`crash-verification.json`、原 manifest hashes、
`measurements.json`、binary 封存與來源／證據 hashes。各次 wall 包含初始化、
I/O、故障測試與數值計算，不是 checkpoint 純 I/O 時間；沒有用它宣稱效能提升。
CUDA allocation 為 N/A，checkpoint memory 額外含 owned candidate 與既有 replicated arrays。

六個成功 MPI 作業共 18 份 rank reports，退出碼皆 0、stderr 皆空。最大 rank wall
與各 rank peak RSS 如下；typed 與 sanitizer save 有執行重疊，不作效能比較：

| 作業 | 最大 wall s | rank 0／1／2 peak RSS bytes |
|---|---:|---|
| release save | 3.01987 | 43,692,032／45,375,488／45,285,376 |
| release resume | 1.46653 | 42,840,064／45,383,680／44,376,064 |
| release recovery | 1.46636 | 43,061,248／44,318,720／45,367,296 |
| typed regression | 8.93034 | 48,328,704／45,580,288／48,500,736 |
| sanitized save | 45.06342 | 195,174,400／207,159,296／200,867,840 |
| sanitized resume | 20.50199 | 193,482,752／202,313,728／198,823,936 |

typed regression 的 exclusive assembly／solver setup／linear solve／明確量測的
communication phases 另存於 `measurements.json`，不重複加總 nested totals。
本批沒有單獨量測 checkpoint 純 I/O latency 或大型網格的記憶體 scaling。

## 剩餘工作

下一步接入 0D／1D／3D native graph providers、完整 configuration／execution identity、
accepted history prefix、pressure guess／donor 與全作業 candidate publication，加入
CLI 保存／恢復與 commit 前中斷驗收。此批的 fixture identity 涵蓋輸入 bytes 與模式，
不是 production graph identity builder。其後才驗收完整 graph 的三個中斷點；VCA、
immersed／moving／FSI、不同 rank 數與跨節點驗收亦未在此批完成。整份 TODO 維持
14／38 已完成，05C 為部分完成。
