# HPC-06A／B／C 完成報告

狀態：完成。日期：2026-09-11。驗收環境：單一 workstation、ParaView 5.13、
PETSc 3.15、Open MPI，最多 4 ranks。跨節點檔案系統與排程器驗收仍屬 HPC-09。

## 完成範圍

HPC-06A 的 source inventory 已盤點完整場 gather、資料庫 required-element 索引、
移動幾何雙 epoch、Bezier mesh 與輸出暫存的大小及生命週期。Flow writer 移除了
一份 4N doubles 的輸出陣列副本；`--memory-report` 會在所有 ranks 收集目前 RSS、
累積 RSS 高水位及 PETSc counters。1024 元素、2299 控制點、4-rank 的 flow
HDF／PVTU 求解都保存各 rank 的 phase records 與 wrapper peak RSS。

HPC-06B 的 PVTU 路徑以 owned elements 建立 cubic Bezier pieces，透過完整
extraction signature 選定共享 Int64 point identity，只交換所需控制 rows 與代表
tuple。Root 只發布 O(ranks) piece filenames 的 PVTU/PVD metadata，不收集完整解場。
CPU flow 與 configured transport CLI 均支援 `--visualization-format pvtu`；既有
VTU／VTKHDF 路徑及 checkpoint 格式保留。初始化的全域 geometry 認證仍在 root，
完成後會釋放該 mesh。

HPC-06C 新增 `--diagnostic-every N`，並保留 `--output-every N` 與
`--checkpoint-every N` 的獨立語意。三者預設值維持既有行為。Flow 與 transport
各以 1／2 ranks、頻率 1／2 驗證；另有 transport 的 PVD 暫存路徑阻擋及新路徑
健康重試。現有本機資料顯示降低頻率會同時降低檔案數與寫入 bytes，未顯示需要
新增有限 I/O aggregators，因此保留每 rank 一份 piece 的設計。

## 數值與資源驗收

1024 元素案例的 HDF／PVTU 各四份 rank reports 都成功，ParaView 讀取三幀；
場最大 relative L2 為 `2.2544463435519334e-15`。HDF 的各 rank peak RSS bytes
為 `348966912, 412311552, 326963200, 308318208`，PVTU 為
`342953984, 410812416, 326205440, 307474432`。三幀 HDF 是 1 檔、1,538,011 bytes；
PVTU 是 16 檔、14,583,075 bytes。PVTU 使用 ASCII 且每幀重複幾何，這組數據不支持
以 PVTU 取代壓縮 HDF 作為單 rank／共享檔案系統的通用預設，也不支持宣稱其降低
整體 peak RSS；其用途是移除大型解場的 root gather 並允許並行 piece 寫入。

頻率矩陣共 9 個成功作業與 1 個預期拒絕。頻率從 1 改為 2 時：

| 路徑 | ranks | 可視化檔案數 | 可視化 bytes |
|---|---:|---:|---:|
| flow | 1 | 7 → 5 | 25,148 → 15,372 |
| flow | 2 | 10 → 7 | 29,553 → 18,307 |
| transport | 1 | 7 → 5 | 19,857 → 12,540 |
| transport | 2 | 10 → 7 | 24,257 → 15,481 |

每組 dense／sparse 的 final checkpoint SHA256 相同。ParaView 讀取 transport 的
10 個時間幀，`three_red`／`three_blue` 在 rank 與頻率配置間的最大 relative L2
為 `1.6766462693301943e-12`，低於 `1e-6` 門檻；shared tuple、GlobalPointIds 與
GlobalCellIds 亦通過。PVD pending path 被目錄阻擋時，兩個 ranks 都在 timeout 前
非零退出、沒有 final PVD 或成功摘要；換新輸出路徑後成功，final checkpoint 與
健康基準相同。

## 重現與證據

```bash
python3 scripts/hpc_flow_memory_regression.py \
  --reference-root outputs/hpc06/flow-pvtu-v1/cli \
  --output-dir NEW_MEMORY_OUTPUT
python3 scripts/hpc_io_control_regression.py \
  --flow-fixture-root outputs/hpc06/flow-pvtu-v1/cli \
  --transport-fixture outputs/hpc01/solver-stdout/native-driven/transport-input \
  --output-dir NEW_IO_OUTPUT
pvpython scripts/test_transport_pvtu_paraview.py NEW_IO_OUTPUT \
  --output NEW_IO_OUTPUT/paraview.json
python3 scripts/hpc_io_scaling_completion.py --output-dir NEW_COMPLETION_OUTPUT
```

最終 completion audit 為 `outputs/hpc06/completion-abc-v2/acceptance.json`，會驗證
12 份 component、flow、transport、large-case 與 reader 證據，並記錄其 SHA256
及目前 production sources。原始大型證據為
`outputs/hpc06/duct-output-large-v1/audit.json`；當前 CLI 頻率與 transport reader
證據為 `outputs/hpc06/io-control-v4/summary.json` 與 `paraview.json`。

這項完成狀態不包含 process-loss recovery、fsync durability、跨節點 metadata
server 壓力或 root 全域 geometry 認證的分散化。這些限制不會讓大型解場重新在
輸出時 gather 到 root；跨節點可攜性與檔案系統行為由 HPC-09 驗收。
