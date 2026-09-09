# HPC-03C：分散式固定幾何 velocity history

日期：2026-09-09。基準 `8f78e2a` 加本批修改。
狀態：完成獨立 history 元件與 backward-Euler 元素接點。HPC-03C 仍部分完成；
完整暫態 operator、Newton、守恆／port 診斷與 graph 尚待接入，保留模式限制。

## 實作

[ImmersedDistributedVelocityHistory.hpp](../../solvers/cpu/include/ImmersedDistributedVelocityHistory.hpp)
使用既有 `ImmersedActiveLayout` 的穩定 node IDs。來源仍為 owned 四場 PETSc
Vec；只抽出每個 owned node 的三個速度，另以 VecScatter 交換本 rank 指定的
節點。壓力、port multiplier 與 gauge 不屬於 history。符號 layout 可複製，
沒有逐 rank 複製完整數值場。來源 ownership 必須以完整四場節點為單位，
可容許零列、零積分工作，以及持有全部列但沒有積分工作的 rank。

`Freeze` 共同檢查 layout、communicator、ownership、有限且前進的精確時間、
相鄰 step index，再抽取 committed velocity 並交換 halo。所有資料 owner 的
速度都必須有限，即使該 rank 沒有積分工作。所有可能失敗的操作完成後才
發布 frozen halo／clock；失敗保持 idle，允許同一物件重試。已凍結期間拒絕
再次 Freeze。Newton 改變來源 Vec 不會改變 history。

`Localize` 依元素 global connectivity 查找 halo，並要求 assembly time 完全
相同；缺少節點、空元素、未凍結與 closed 狀態都拒絕。新增具名 transient
元素 overload 使用相同的 backward-Euler kernel，並核對 dt 對應 target time。
既有序列 overload 與資料格式不變。

`ReleaseTrial` 是 runtime 在共同 commit／abort 後於每個 rank 執行的本地
noexcept 操作；它只取消 snapshot，不會自行提交流場或推進 accepted clock。
Close 可重複呼叫，會清除 active 狀態並釋放所有 PETSc handles。

數值 buffer 包含 3×owned nodes 的分散 Vec，以及 3×required nodes 的 halo
Vec、candidate 與 frozen 陣列。兩份陣列用於失敗前不發布的交換，不能只用
單一 frozen 陣列的大小代表元件記憶體。沒有宣稱本批已完成 force freezing、
完整暫態求解、moving history extension 或跨節點驗收。

## 驗證

```bash
make -C solvers/cpu immersed_distributed_history_test CXX=mpicxx \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
python3 scripts/hpc_immersed_history_regression.py \
  --output-dir outputs/hpc03/distributed-history/final --split
```

硬體／軟體：本機 TsungYehLab、GCC 11.4、OpenMPI 4.1.2、PETSc 3.15.5
real64/int32；OMP／BLAS threads=1。MPI 使用允許本機 socket 的執行權限。
這是小型本機 correctness test，不是 scheduler 大案例或效能擴展驗收。

測試為完整 unit cube、4 cells、112 active nodes、兩個 controller 加 gauge。
每個 communicator 跑一般 ownership、所有工作／列集中 rank 0，以及工作在
rank 0／全部列在最後 rank 三種配置。1／2／4 ranks 加 world 3 的 split 1+2，
共 4 MPI jobs、10 rank reports、30 組 completion observations。

檢查兩份非零 committed velocity 的凍結與精確 clock、逐節點序列 identity
history 比較、backward-Euler conservative Jacobian／residual 比較、相對零
history 的非零殘差差異、來源 Vec 修改隔離、單 rank 錯誤 layout／node ID／
clock／index、不同但各自有效的 dt、active／released／closed 拒絕，以及
只有一個資料 owner 的 NaN 故障後原物件重試。測試刻意使壓力／scalar 為 NaN
來確認 history 僅抽取速度；這不是允許完整 runtime 接受非有限流場。

最終結果與 binary SHA／各 rank wall time、peak RSS 保存在
`outputs/hpc03/distributed-history/final/summary.json`。history 逐節點與元素比較
要求完全相同，並要求 history 對殘差的最大影響 > 1e-6。

初次 fixture 使用 compact catalog 呼叫只接受展開 quadrature 的序列參考 API，
因此在數值比較前拒絕；已修為 expanded fixture，原紀錄保留於 `initial/`。
`expanded/` 是加入最後單一 owner 故障案例前的通過結果；以 `final/` 驗收。

最終四個作業均 exit 0、無 timeout。30 組 observation 的元素最大差為 0，
history 對殘差最大影響為 0.00242323。最大 rank wall time
為 1.945833 s，peak RSS 為
36818944 bytes；這些包含 fixture、參考計算及啟動，並非線性求解時間。
一般分區的 2／4 ranks 各自 owned nodes 與 halo nodes 均小於 112。
