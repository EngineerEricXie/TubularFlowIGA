# HPC-03C：owned 暫態體積組裝與凍結外力

日期：2026-09-09。基準 `bf29fd1` 加本批修改。
狀態：完成固定幾何 backward-Euler 體積組裝元件。完整暫態 operator、Newton、
port／wall／守恆診斷及 graph 接入仍待接續，HPC-03C 保持未勾選。

## 實作

[ImmersedDistributedTransientVolume.hpp](../../solvers/cpu/include/ImmersedDistributedTransientVolume.hpp)
接受 immutable catalogs、`ImmersedActiveLayout`、owned state 模板及本 rank
owned cells。共同驗證 geometry／layout／quadrature 設定及每個 active cell
恰有一個積分 owner；缺少、重複或無效 ownership 在建立 history halo 前拒絕。
owner 與 state rows 可以不同，包含完全空 rank 與純遠端資料積分。

`Freeze` 先共同驗證 rho／mu／dt 與 evaluator，再凍結上一時間步 velocity。
只在 cell owner 依 canonical volume-point 順序評估 body force，每個點恰一次。
任一 callback 拋例外或回傳非有限數值時，共同撤回新 history 和 force candidate；
原物件可重試。全部成功後才發布外力與 transient parameters。外力 callback
必須是本地操作，並在各 rank 代表同一物理場；函式本體不能用文字一致性檢查。

`BuildVolume` 是可放入 `ImmersedDistributedAssembly` callback 的本地操作。
它依 global connectivity 取得凍結舊速度，用對應 cell 的凍結外力重播既有
conservative backward-Euler kernel，並檢查點數完整消耗。後續 callback 或
來源 Vec 的修改不會改變輸入。外力僅儲存 owned cells 的三分量點值，不複製
全域數值場。compact／expanded 使用各自既有 canonical point iterator。

`ReleaseTrial` 是共同決定後各 rank 執行的本地 noexcept 撤回操作。Close 清空
active snapshot 並釋放 history handles，可重複呼叫。本元件不自行提交解、
推進 accepted clock、認證 material wall velocity，也不建造 wall／port／ghost／
gauge 項。這些仍須在完整暫態 operator／runtime 接入時一起驗證。

## 驗證

```bash
make -C solvers/cpu immersed_distributed_transient_volume_test CXX=mpicxx \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
python3 scripts/hpc_immersed_transient_volume_regression.py \
  --output-dir outputs/hpc03/transient-volume/final --split
```

本機 TsungYehLab；GCC 11.4、OpenMPI 4.1.2、PETSc 3.15.5 real64/int32，
OMP／BLAS threads=1，MPI 在允許 socket 的本機環境執行。

案例是背景 `[0,1]^3` 的 4×1×1 cells、實體 `[0.125,0.875]^3` 切割立方體，
112 active nodes、451 rows。scalar rows 在本 volume-only 測試保留結構零，
沒有解這個缺少完整邊界項的系統，也沒有把 exit 0 解讀為 Newton 收斂。

1／2／4 ranks 與 world 3 的 split 1+2，各跑一般分區、全部工作／列在 rank 0、
工作在 rank 0 而全部列在最後 rank 三種配置，各搭配 compact／expanded。
每次將本地 volume block 插入真正的 distributed Mat／Vec，再比較所有 owned
negative residual 與 Jacobian action，參考是既有序列 global-ID history 及
expanded kernel 的全域積分。門檻為最大絕對差 < 1e-12；非零舊速度相對零
history 的殘差差異必須 > 1e-6。

另檢查兩份非零 history／外力、缺少與重複 owner、不同 catalog mode／mu、
callback 拋例外／NaN、凍結後禁止再次執行 callback、部分插入後組裝失敗與
原物件重試、空 rank 無外力資料、Close 後拒絕。失敗重試的相同輸入重用
已算好的序列基準，仍逐次比較全部 owned 殘差與 action；新輸入重新計算基準。

每個 rank 的 `integration_s`、`halo_s`、`stash_s` 來自最後成功的 assembly。
`run.json` 的 wall time／peak RSS 包含 fixture 與序列參考成本，不代表 solve
時間或擴展效能。沒有線性求解階段，也沒有 GPU 或跨節點測試。

最終驗收保存在 `outputs/hpc03/transient-volume/final/summary.json`，包含 binary
SHA、MPI 命令、每個 rank 的退出狀態、資源報告與數值 observation。
`initial/` 為補上非有限外力及單 owner 部分插入案例前的第一版測試結果。

最終 4 MPI jobs、10 rank reports、60 observations 全數通過，均 exit 0、
無 timeout。最大 residual／Jacobian-action 絕對差 4.5976e-13；非零 history
對殘差的最大影響 0.00201951。最大 rank wall time 100.301733 s，
peak RSS 50012160 bytes。編譯無新增 compiler warning；受限編譯環境的
`opal_ifinit socket errno=1` 不影響這些已通過的授權 MPI runtime 測試。
