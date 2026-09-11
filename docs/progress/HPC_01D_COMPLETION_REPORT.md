# HPC-01D：執行資源與能力檢查完成報告

日期：2026-09-11。基準 revision `b74b1ccb08ef06690ae3ddb54c0f0fc99f84e16c`
加本批兩個入口補強。

## 完成範圍

所有 production PETSc 路徑都在讀取大型輸入、建立 solver 或輸出前進入資源
preflight：CPU flow、configured／legacy transport、mesh check、assembly smoke、
native 1D、multidomain／bifurcation、sequential 與單 rank FSI exporter。最後兩個
缺口已在本批接入 `RequireExecutionResources`。CUDA 五個入口使用對等的單程序
`RequireExecutionEnvironment`，並在 GPU 查詢前拒絕不支援的 launcher／thread 設定。

共同 CPU 檢查會：

- 驗證 communicator rank 數與資料庫 partition 數，並以 `PetscInt` 上限檢查
  node×field rows；
- 拒絕非 real64 PETSc，逐 rank 比較 PetscInt／PetscScalar ABI，摘要明列 index
  與 scalar 位寬；容量公式直接使用該 build 的 `numeric_limits<PetscInt>`；
- 解析 OpenMP、thread limit 與 OpenBLAS／MKL／BLIS requests；多 OpenMP threads
  要求實際 MPI thread level 至少 FUNNELED；
- 顯示 MUMPS／Hypre build flags，並在明確選擇 LU／Cholesky backend 時，以實際
  matrix 查詢可用性及比較全 rank 的有效 PC／backend。

未指定 factor backend 時維持 PETSc 的選擇權；KSP／PC 的完整選項與 sub-KSP
調校屬 HPC-04。HDF5 writer 的 schema、close 與容量屬 HPC-06；scheduler allocation、
CPU／NUMA／GPU binding 與跨節點證據屬 HPC-09。這些分工不構成 HPC-01D 缺口。

## 最終驗收

環境為 GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real64／Int32／MUMPS，
OMP／OpenBLAS 各 1 thread，除異質 thread 案例外。

| 驗收 | 結果 |
|---|---|
| Resource unit | world3 與 split1+2 通過；17／16 個案例，含實際 MPI thread level、ABI agreement、容量、report stream 與 factor queries |
| 三種 CPU MPI 工具 | 39 個原生作業及 geometry world／split groups 通過；partition、field capacity、backend 與失敗摘要均受檢查 |
| CUDA host／native／數值 | 84 個 host checks、23 個原生啟動、38 個數值／stdout 作業與 27 個場比較通過；實際 2-rank Open MPI launcher 被拒絕 |
| 最終 production entry matrix | 20 個作業、六個 binaries、九個輸入 hash 通過；12 個 thread 拒絕、四個 database 拒絕、backend 拒絕及三個健康作業符合預期 |
| 新版 coupling 回歸 | 80 個 graph／bifurcation／sequential 作業及 222 份輸出比較通過；28 個 cleanup failure 未發布成功結果 |

最終入口矩陣在 `outputs/hpc01/resources/final-entry-matrix/summary.json`；coupling
回歸在 `outputs/hpc01/resources/final-coupling-regression/summary.json`。兩個新增
入口與 resource unit 均強制重建，沒有 compiler warning；`git diff --check` 通過。

## 能力界線

本機只有 PETSc real64／Int32，因此沒有宣稱實跑 64-bit-index 或 complex build。
Complex build 會由明確 ABI gate 拒絕；64-bit index build 會在摘要顯示 64 並以其
實際上限作容量檢查。跨 build 的編譯與數值矩陣由 HPC-09 的平台建置驗收追蹤。
資源摘要描述 runtime 可見設定與編譯能力，不等同實測 worker 數、核心綁定、
效能或跨節點擴展性。
