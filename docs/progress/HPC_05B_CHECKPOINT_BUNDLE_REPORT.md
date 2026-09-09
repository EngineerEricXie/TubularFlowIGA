# HPC-05B：版本化 checkpoint bundle 與發布驗收

日期：2026-09-09。基準 revision：`cd03c1818698336da9b9be19708b7ab023ea14de`
加本批修改。狀態：**HPC-05B 已完成**，範圍為獨立 bundle 格式、發布、驗證與
恢復搜尋。HPC-05C 的 domain capture／restore 和 MPI graph 續跑尚未完成。

## 實作範圍

- [Manifest](../../include/CoupledCheckpointManifest.hpp)：版本 1、完整 metadata
  checksum、精確 uint64／binary64、大小／欄位／catalog／owner／配置嚴格檢查。
- [Bundle](../../include/CoupledCheckpointBundle.hpp)：immutable generation、暫存
  shard、writer receipt、落盤資料重新校驗、file／directory fsync、完成 manifest
  最後發布；串流載入與最後完整且相容 generation 的搜尋。
- [協議文件](../architecture/COUPLED_CHECKPOINT_BUNDLE.md)：API 順序、caller 的
  accepted-boundary／完整 catalog／MPI agreement 責任、失敗與恢復語義。
- [測試](../../solvers/coupling/tests/test_coupled_checkpoint_bundle.cpp)：真正的
  fork producer／新 loader 程序、程序終止、資料破壞、故障注入與 32 MiB 串流。

舊 `.ntiga`、flow／transport／1D／VCA checkpoint 格式與 numerical runtime 沒有改動。
本批新增的 payload 容器不把 visualization snapshot 或舊 standalone state 自動轉成
可恢復的 graph。執行模式為單工作站 POSIX shared filesystem、兩個並行 producer，
由 parent 收齊 pipe receipts 後發布；沒有呼叫 MPI 或 CUDA。

## 驗收標準與結果

預先要求完整寫入可由獨立程序讀回，所有損毀／不相容案例拒絕；中斷發布不覆寫舊版本；
不同 epoch 的分片不能通過；parser 有明確上限且精確往返 clock／uint64；串流記憶體
不隨 payload 大小配置一份完整場。以下均由測試中的實際判定驗證：

| 案例 | 結果 |
|---|---|
| 兩個並行 producer、binary payload 含 NUL／newline、新程序載入 | 全部 bytes 與 catalog 一致 |
| manifest 每一個截斷位置、尾端雜訊、錯版本／負數／整數溢位／非 canonical 數字、超限 | 全部拒絕 |
| uint64 最大值、0.1 累加的實際 binary64 clock | 精確往返，不以 `N*dt` 重建 |
| NaN／infinity／零／負時間、重複／非法 ID、owner 超出 communicator | 全部拒絕 |
| 寫分片中 `_exit(87)`、manifest 暫存已 sync 後 `_exit(86)` | 無完成 marker；新的 loader 程序找回前一完整 epoch |
| writer receipt 後損毀分片、缺 receipt、寫入／fsync 注入失敗、producer exception、寫多／寫少 | 發布失敗，舊 manifest 位元組保持一致 |
| 完成後 byte corruption／截斷／多餘 bytes／缺片／symlink／混 epoch 分片 | 明確載入拒絕；搜尋回復前一完整相容 epoch |
| case／configuration／execution／rank／format 不相容 | 全部拒絕 |
| 完成 manifest 損毀、同最大 step 的兩個完成 generation | 前者拒絕；後者報歧義，仍可明確選擇有效 epoch |
| 空 rank 的 zero-byte shard、32 MiB 分片 | 都可發布、載入與串流驗證；大分片逐 byte 比較 |

最終 optimized 與 sanitizer 各通過 **734 個預期拒絕**；總數包含每個 manifest
截斷位置，不代表 734 次 MPI 作業。程序終止的預期退出碼為 86、87，其餘 producer／
loader 與整體測試退出碼為 0。

## 命令與環境

Linux x86_64 工作站、GCC 11.4、C++17，`-Wall -Wextra -Wpedantic` 無新警告。
使用新的、尚不存在的 output child directory，每次重跑更換末段名稱：

```bash
make -C solvers/coupling coupled_checkpoint_bundle_test CXX=g++
mkdir -p outputs/hpc05/bundle
solvers/coupling/coupled_checkpoint_bundle_test "$PWD/outputs/hpc05/bundle/final"
g++ -Iinclude -Isolvers/cpu/include -std=c++17 -O1 -g \
  -Wall -Wextra -Wpedantic -DIGA_COUPLED_CHECKPOINT_TESTING \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  solvers/coupling/tests/test_coupled_checkpoint_bundle.cpp \
  -o outputs/hpc05/bundle/coupled_checkpoint_sanitizer_test
ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=halt_on_error=1 \
  outputs/hpc05/bundle/coupled_checkpoint_sanitizer_test \
  "$PWD/outputs/hpc05/bundle/sanitized-final"
git diff --check
```

亦有 opt-in `make -C solvers/coupling checkpoint-bundle-test`，在新 `/tmp` 子目錄執行。
測試用故障 hooks 只在 `IGA_COUPLED_CHECKPOINT_TESTING` 定義時編譯。

權威證據位於 `outputs/hpc05/bundle/`：`build.log`、`final.log`、
`sanitized-final.log`、空的 `sanitized-final.err`、封存 test executables 與
`acceptance.json` 的來源／執行檔 SHA-256。測試生成的 immutable／partial／corrupt
bundles 保留在對應 output 目錄供檢查，不提交到 Git。

最初 restricted sanitizer run 的 LeakSanitizer 因 sandbox ptrace 限制退出 1；
取得所需執行權限後，最終版 ASan／UBSan／LeakSanitizer 全部通過、退出 0、stderr 0 bytes。
`run1`／`run2`／`sanitized` 為較早版本的中間檢查，不作為最終來源的替代證據。

## 時間與記憶體

optimized 最終版，32 MiB payload 加一個 empty shard：

| 項目 | 量測 |
|---|---:|
| streaming write＋file／directory sync | 0.130980 s |
| 依 writer receipt 重新校驗＋manifest write／publish／sync | 0.0899502 s |
| 完整 load verification＋第二次 streaming read／逐 byte 判定 | 0.180584 s |
| 全部格式／破壞／中斷／I/O 測試 wall time | 0.510068 s |
| parent peak RSS | 8000 KiB |
| child 程序中最大的 peak RSS | 3180 KiB |

以 `steady_clock` 與 Linux `getrusage(RUSAGE_SELF/RUSAGE_CHILDREN)` 量測。
sanitizer 的額外記憶體成本不作效能基準。這是本機、含 cache 影響的小型 I/O 證據，
沒有宣稱 shared parallel filesystem、跨節點吞吐或 rank scaling。
assembly／solve、numerical relative L2、convergence、守恆與 CUDA allocation 為 N/A：
本項測試的是無 numerical solve 的格式及發布元件，場解／全歷史 restart 由 05C/D 驗收。

## 剩餘工作

HPC-05C 應從 validated graph 建立完整 catalog 與 compatibility digests，在 accepted
macro-step 邊界擷取 runtime／pressure／donor／history，跨 MPI ranks 收齊 receipts
及錯誤，再呼叫本協議。載入時先構造 unpublished candidates，全部 domain 驗證後才
一起恢復 accepted state。仍需 0D／1D／貼體 3D 的不中斷／新 MPI job 續跑全歷史比較，
含真正的 graph commit 前中斷。05D 加入 immersed／moving／FSI 與不同 rank ownership。

coordinator 目前以 64 KiB buffer 串行驗證每片；06 可量測後改善分散式校驗／I/O。
fsync/link 的 power-loss 與跨節點儲存保證尚未作硬體故障測試，不能從本批程序終止
驗收推導。舊 generation retention、排程 signal 的安全 step-boundary 綁定亦待整合。
