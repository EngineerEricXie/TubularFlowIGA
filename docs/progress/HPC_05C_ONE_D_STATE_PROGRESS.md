# HPC-05C：1D accepted state 與串流 checkpoint 進度

日期：2026-09-09。基準 revision：`9d0676283e1dfd2c6fd159651ac0c505998304e2`
加本批修改。狀態：**部分完成，HPC-05C 保持未勾選**。

本批完成原生 1D runtime 的 accepted capture／fresh restore，以及可與 05B bundle
搭配的 metadata、field 兩分片格式。五種模式從第 3 個 macro-step 的磁碟 bundle
恢復，新的 executable 接續後面 7 步，與父程序預先寫下的不中斷參考逐位元一致。
正式 native graph CLI、1D provider wiring、3D providers 與完整歷史 prefix 尚待實作。

## 狀態介面與格式

[OneDRuntime.hpp](../../solvers/one_d/include/OneDRuntime.hpp) 新增
`OneDFlowCheckpointState`、`CaptureCheckpointState()`、`RestoreCheckpointState()`：

- 擷取完整 hydraulic arrays／outlets、species arrays／boundary flux／step accounting、
  `LastInlet`、動態 segment radius、四個 physiology／perfusate Hct／Hb scalar。
  保存 accepted macro count、最後 start／dt，與既有 configured step count、physical
  time、64-bit internal substep count。計數只在 finalize 增加，abort 不覆寫 accepted clock。
- 僅接受固定 macro dt、整數 configured subcycling。這與目前 graph 契約一致；
  variable macro dt、multirate 不在本批支援範圍。先檢查 dt ratio 範圍再 `llround`，
  防止不可信 checkpoint 的巨大 dt 進入不合法的整數轉換。
- constructor 可選參數為 provider 驗證過的配置 identity SHA-256，必須涵蓋完整配置、
  selected flow system、network、外部輸入與其順序。此介面不會替 caller 建立 digest。
  舊 caller 可以省略 identity，但使用新 checkpoint API 時會明確拒絕。
- capture 要求 Ready 且已有 accepted step。restore 僅限 fresh、initialized Ready
  candidate。先驗證身分、形狀、clock、outlet／species models、accounting coverage、
  有限值，再重建 area0／resistance；最後以不丟例外的 moves／swap 發布。
  驗證失敗不修改 fresh runtime；已運行的 owner 禁止直接覆寫。
- 沿用原傳輸 solver 的數值語義，有限的負 interior concentration／initial mass 可保存。
  現有正值 area／radius 與非負 inlet concentration 要求保持不變。

本批修正 05A 盤點中的遺漏：`ApplyOneDVasodilation` 使用 const configuration，
但 `ApplyOneDCoupledInlet` 會改寫 physiology 和 coupling perfusate oxygen 的 Hct／Hb；
也會改寫 `inlet_value` 並清除該 species 的 waveform。後續 inlet 可以省略這些欄位，
所以僅保存 `LastInlet()` 不足以恢復。現在會保存這些動態值；配置中的原始 waveform
名稱另保留作 restore 驗證，不因 `InitializeCoupled` 清除當前 waveform 而遺失。

[OneDAcceptedCheckpoint.hpp](../../solvers/one_d/include/OneDAcceptedCheckpoint.hpp)
定義 `IGA_ONE_D_ACCEPTED/1`：

| 分片 | 內容與限制 |
|---|---|
| metadata | identity、clock、64-bit counters、完整 inlet presence／species／metadata、blood state、outlet 動態值、species inlet／flux／accounting、每個 field 長度；使用共用 16 MiB metadata 上限及 4096-byte strings |
| fields | 依序為 area、flow、pressure、node pressure、segment flow、segment radii、每個 species concentration；IEEE binary64 little endian，最多 64 KiB 每次 sink 呼叫 |

解析 metadata 時，storage shape 只來自身分已核對的初始化 target，不依不可信長度
配置大型陣列。串流 consumer 接受任意 byte fragmentation，只有 8-byte partial scalar
與 vector pointers 等少量額外 storage；拒絕多餘、截斷與非有限場值，失敗後不可重用。
metadata 的 active inlet scalar 必須有限，inactive storage 的原始 bits 仍完整保留，
包含 signed zero／NaN；inlet metadata 的空字串 key 也按原有介面保留。

兩片必須同屬一個完整 bundle。caller 必須在所有分片 SHA／size、consumer `Finish()`、
runtime validation 與群組 agreement 通過後，才把候選 owners 換入 live graph。
codec 本身沒有 MPI；新的 test consumer 只修改尚未發布的 snapshot。現有 legacy
`OneDCheckpoint` 及 `RestoreCommittedState` 格式／呼叫行為保持相容。

## 案例與驗收結果

[測試程式](../../solvers/one_d/tests/test_one_d_accepted_checkpoint.cpp) 產生兩條 root
branches、兩個 outlet、每段兩個 cells 的 SWC fixture，運輸 signal。configured dt=0.001、
macro dt=0.002、horizon=20 configured steps。配置／network／waveform 的內容與測試
mode recipe 綁定 identity。比較門檻預先設為新格式逐位元一致，包括 derived geometry、
完整 metadata／fields、三個 port 與 species accounting。

| 驗收 | 結果 |
|---|---|
| coupled vasodilation | 非零 radius 變化、Hct=35／Hb=12.25；保存時 LastInlet 已省略 Hct／species，但前步 blood state 與清除的 waveform 仍精確恢復 |
| staged transport | hydraulic frames 在下一步重建，既有 scalar accounting 與後續正／反向流一致；不啟用該路徑原本拒絕的 vasodilation |
| RCR | reference／capacitor pressure、outlet flow 與後續動態狀態一致 |
| compliant explicit | area／flow／pressure 及 internal substeps 精確恢復，後續 explicit 更新一致 |
| configured waveform | 保留未清除的 waveform，從正確 clock 接續評估外部 periodic table |
| 新程序 | 五個 mode 各讀取磁碟兩分片 bundle，與父程序留下的第 3–10 步八份完整 state fingerprints 一致，共 35 個續跑 macro-steps |
| 負面／lifecycle | 每個 metadata 與 field 截斷位置、尾端 bytes、非有限場、錯 identity／clock／catalog／model／coverage、非法半徑、active capture、重複 restore、不同 dt trial abort、prepared abort 均驗證 |
| 寬計數與 signed state | 超過 `INT_MAX` 的 internal substeps 可 roundtrip；有限 signed interior state 不被新增 positivity gate 改寫 |
| 大型 payload | 24 MiB field roundtrip 通過，單次 chunk 最大 65536 bytes，無 16 MiB field 限制；只證明 streaming 元件，未宣稱大型 graph scalability |
| 既有純 C++ tests | core、coupling、runtime 均通過 |
| MPI failure regressions | trial 與 adapter 各一個三 rank job；WORLD 與獨立 1／2-rank groups 的 local failure、abort／retry、COMM_SELF 參考比較均通過 |
| legacy checkpoint | 25 個三 rank jobs（75 份 rank reports）通過，含 rigid／explicit／species 續跑、PETSc binary bytes 相容性，以及 14 類讀寫故障的一致拒絕；數值場最大誤差 0，clock 的最大 relative L2 約 6.782e-17 |

legacy gate 沿用 `hpc_one_d_cli_regression.compare`：非零參考 norm 的 relative L2
不超過 1e-6，零參考的絕對 error 不超過 1e-12；restart clock 逐列使用 1e-12
相對尺度門檻。ASan／UBSan／
LeakSanitizer 使用同一新測試 source，包含 fresh executable 與大型 stream。
optimized 與 sanitizer 最終測試各通過 4,573 個預期拒絕判定，退出碼均為 0；
ASan／UBSan／LeakSanitizer 無錯誤診斷。sanitizer peak RSS 為 121,768 KiB。

## 環境、命令與證據

Linux x86_64 工作站、GCC 11.4、OpenMPI 4.1.2、PETSc 3.15.5 real64/int32，
1 thread；MPI failure tests 使用 preonly／LU／MUMPS，未固定 CPU affinity。
編譯開啟 C++17、`-Wall -Wextra -Wpedantic`，無新警告。小型本機驗證不涉及 cluster
login node 或 GPU；本批不宣稱 MPI checkpoint 整合／跨節點／scaling。

```bash
make -C solvers/one_d one_d_accepted_checkpoint_test core-test CXX=g++
mkdir -p outputs/hpc05/one-d-state
OMP_NUM_THREADS=1 /usr/bin/time -v solvers/one_d/one_d_accepted_checkpoint_test \
  outputs/hpc05/one-d-state/final-v2
make -C solvers/one_d one_d_trial_failure_test one_d_adapter_failure_test \
  iga_1d one_d_checkpoint_format_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
# 為 trial／adapter 分別建立新的 --output-dir，再執行：
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
PETSC_OPTIONS='-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps' \
mpiexec --oversubscribe -np 3 python3 scripts/hpc_rank_run.py \
  --output-dir outputs/hpc05/one-d-state/trial-mpi --expected-ranks 3 --timeout 120 \
  -- solvers/one_d/one_d_trial_failure_test examples/one_d/multispecies_physiology
python3 scripts/hpc_one_d_checkpoint_regression.py \
  --output-dir outputs/hpc05/one-d-state/legacy-regression
```

每次重跑必須使用新的 child directory。sanitizer flags 為 `-O1 -g
-fsanitize=address,undefined -fno-omit-frame-pointer`；執行設定
`ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=halt_on_error=1`。MPI 與 LeakSanitizer
均使用所需本機執行權限。restricted build 的 `opal_ifinit errno=1` 屬環境訊息，
實際 MPI jobs 皆成功；沒有將 sandbox 限制當成缺少硬體。

證據位於 `outputs/hpc05/one-d-state/`：最終 `final-v2.log/.err`、
`sanitized-v2.log/.err`、`trial-mpi`／`adapter-mpi` rank reports、legacy regression
summary、`acceptance.json` 的 85 份來源與 1,045 份證據／輸入／binary 雜湊、
封存 executable，以及 commit 後的 `commit.json`。前期 fixture 缺 reverse outlet
concentration、缺 waveform definition 的兩次失敗保留為 attempt-1／2；已修正後通過，
不把失敗納入成功統計。

新測試分列各 mode 的 native solve 呼叫總時間、bundle write／publish、
load／verify／restore 與 24 MiB stream encode／decode 時間。這是單次正確性測試的
計時，沒有 speedup 主張；native 小型 1D solve 未使用 PETSc matrix assembly，
assembly／preconditioner timing 為 N/A。MPI failure tests 的 wall 含啟動／故障注入，
不可當作 assembly／solve performance。

optimized 元件測試總 wall 0.14 s、peak RSS 54,232 KiB；包含兩個 field 陣列的
24 MiB stream encode／decode 為 0.02512 s。五個小型 bundle 的 durable write／publish
各約 0.01561–0.01749 s，load／verify／restore 各約 0.000162–0.000259 s。
各 mode 恢復後七次 native solve 累計 0.0000164–0.00008292 s；細項保留於日誌。

trial MPI 的每 rank peak RSS 為 39,141,376／39,800,832／40,509,440 bytes，wall 最大
0.41525 s；adapter 為 40,955,904／40,960,000／40,964,096 bytes，wall 最大 0.51447 s。
這些不是各 rank 同時的 aggregate peak。新元件測試另外配置 24 MiB source 與
24 MiB target，因此 process RSS 包含兩個真實陣列，不能把它誤稱為 codec buffer。

## 剩餘工作與下一步

HPC-05C 仍需要：為 graph 建立完整 compatibility identity 與 1D provider wiring；
補上 body-fitted 3D flow／transport accepted state、provider 與 restored transport count；
保存所有既有 accepted CSV/history prefix；native runner 接上相同 rank bundle、
全作業候選恢復 agreement、續跑 clock、pressure controls 與 donor；最後驗證完整
0D／1D／3D graph 新 MPI job 的歷史及 commit／write／manifest 中斷恢復。

本批沒有把 native graph 的 checkpoint CLI 放行，也沒有完成 moving／immersed／FSI
state providers 或 rank 重分區。完整 TODO 仍為 14／38 已完成，05C 保持部分完成。
