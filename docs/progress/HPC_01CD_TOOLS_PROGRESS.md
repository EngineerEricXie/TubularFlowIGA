# HPC-01C／D：其餘 CPU MPI 工具

日期：2026-09-08。狀態：本批通過；HPC-01C／D 與完整清單仍未完成。
接續 [資源檢查報告](HPC_01D_PROGRESS.md)。

## 變更

`iga_mesh_check`、`iga_assembly_smoke` 與 legacy `iga_transport` 接入共用的
execution resource preflight，保留原有 positional CLI。Database 的開啟、
partition／rank／row capacity 檢查採 collective local stage；單 rank 缺檔、
參數錯誤或分區不合會共同拒絕。Mesh checker 在協調區段內載入 owned elements，
避免其他 rank 先進入幾何 reductions。

Assembly smoke 的 FIELDS 改為完整字串解析，在轉成 PetscInt 前檢查正值與
容量，拒絕 `2junk`、32-bit 環境中的 `4294967297` 等先前可能截斷／誤接受
的輸入。各 rank 的 fields 先精確比較，再建立矩陣；局部元素插入亦加入協調，
assembly、diagonal、matrix info 回傳碼受檢查，已建立矩陣由 RAII owner 清理。

Legacy transport 的案例／參數／標籤／速度／邊界讀取加入共同錯誤處理，
STEPS 使用完整非負整數解析，步數及是否進入 output branch 必須在 ranks 間
一致。保留既有 case parser、物理模型、矩陣與時間步公式。
輸出仍 gather-to-root；scatter 回傳碼受檢查，root read view 與 scatter/vector
使用 RAII，實際檔案寫入 close 後確認成功，再共同報告錯誤。

Legacy KSP 的明確 factor backend 在 setup 前查詢支援能力；setup／solve 回傳碼
及收斂原因加入協調。僅在已檢查回傳碼的 setup／solve 呼叫期間安裝 PETSc
return-error handler，結束後恢復原 handler。沒有更改 KSP defaults 或非零
initial-guess 規則；兩步 KSPPREONLY 因該規則不相容而失敗時，現在共同退出 1。

[InspectGeometry](../../solvers/cpu/include/TransportElement.hpp) 的局部元素／
ownership predicate／quadrature 評估加入 CollectiveLocalStage，再執行原 reductions。
沒有更改積分點、Jacobian 公式、非正 determinant 判定或 bad element/sample 計數。

## 驗收

[工具 regression](../../scripts/hpc_tools_regression.py) 完成 39 項：

| 類別 | 數量與結果 |
|---|---|
| 修改前／後，三工具各 1／2 ranks | 12 項退出 0；Jacobian／matrix 摘要精確相同，legacy 場檔同 rank 逐位元相同 |
| 每個工具的單 rank 缺檔、partition mismatch、參數缺失、非法 threads | 12 項共同拒絕，退出 1 |
| Assembly fields 非法值／溢位／跨 rank 差異 | 6 項退出 1 |
| Legacy 案例／STEPS／output branch 差異、root 寫入及 KSP 錯誤 | 7 項退出 1 |
| 修改前／後的退化 Jacobian | 2 項退出 2；minimum_detJ=0、bad_elements=1、bad_samples=64 維持相同 |

正常 cube 的 minimum_detJ 保持 `0.125`；assembly 沒有 missing diagonal 或額外
matrix mallocs。Legacy 的兩個 fields、64 個節點皆為有限值，整體 norm 非零；1／2 ranks
相對 L2 誤差為 `3.9648210454771405e-13`，通過既有 `1e-6` 門檻。
故障作業沒有輸出成功 summary。各作業 90 秒 timeout，預期退出碼均核對。

[幾何 unit](../../solvers/cpu/tests/test_geometry_inspection_failure.cpp) 在 world3
與獨立 1+2 groups 上驗證：最後一個 rank 的 ownership predicate 丟例外、共同
diagnostic、恢復後重試、單 rank 扁平幾何，以及只有 rank 0 有元素的空 rank
reduction。三群組全部通過，沒有變更幾何數值容許值。

同版 staged／trial／port 的 104／93／64 個案例、flow／species graph 群組數值
比較均通過。VCA runtime 原有預設，以及 MUMPS smoke 的非零流場、reservoir、
1／2 ranks 與同 rank 續跑精確比較通過。六項 MPI runtime 與兩項 VCA 作業各自
在 180 秒時限內退出 0。`make mesh-test` 通過，末尾為 `mesh_core_test: PASS`。

## 基準與修正紀錄

修改前先重建並保存三個工具 binary、原始 sources 與 hashes。
最初從 staged fixture 複製的資料含 boundary label 3，舊版 legacy conversion
只定義 0／1／2，因此 legacy baseline 也拒絕。正式測試在副本上將出口 label 3
映射到 2，同步修改 element face labels 與 VTK point labels；幾何座標不變，
修改前／後 binary 使用完全相同的適配資料。原始 fixture 與失敗 log 保留。

準備腳本起初將 fixture 的格式版本判成 v3，在執行求解器前觸發 assertion；
修正為實際 v5 後重跑。此修正只屬測試資料識別，production `.ntiga` 格式未改。

第一輪完整工具回歸在 KSPPREONLY 故障案例中發現 PETSc 3.15.5 的預設 traceback
handler 直接 abort，沒有返回已檢查的呼叫端。保存該 binary／source／log 後，
加入限定範圍的 return-error handler；以最終 binary 重跑全部 39 項通過。
先前 abort 不算預期退出 1 的驗收證據。

## 重現與證據

環境：WSL 工作站、GCC 11.4、Open MPI 4.1.2、PETSc 3.15.5 real/double、
32-bit PetscInt、MUMPS，OMP／OpenBLAS 各 1 thread。Legacy 正向比較使用
GMRES＋MUMPS，維持既有 initial-guess 行為。
HEAD `ee8da2ab528849814b6d12c8186476b5f7f59124` 加未提交工作樹；建置無新 compiler warning。

```bash
make -C solvers/cpu iga_mesh_check iga_assembly_smoke iga_transport \
  geometry_inspection_failure_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 timeout 180 mpiexec --oversubscribe -np 3 \
  solvers/cpu/geometry_inspection_failure_test
make mesh-test
```

工具比較腳本需要保存的修改前 binary 及包含 `serial.ntiga`、`group.ntiga`、
`controlmesh.vtk`、`initial_velocityfield.txt`、`simulation_parameter.txt` 的
單元素 v5 fixture。這是明確限定的 regression fixture，不是通用 database converter。

```bash
python3 scripts/hpc_tools_regression.py \
  --fixture-source outputs/hpc01/tools/fixture \
  --baseline-dir outputs/hpc01/tools/before \
  --output-dir /tmp/tools-fresh-output
```

Ignored evidence 位於 `outputs/hpc01/tools/`：`before/` 與 `fixture/`、
`preparation-failure.json`／準備腳本快照、`verified-final/` 的 handler 問題紀錄、
`transport-abort-attempt` binary／source，以及最終 `accepted/summary.json`／場檔。
`runtime/mpi-summary.json`、`vca/summary.json`、`build/`、`evidence-summary.json`
與 `inventory.json` 保存命令、退出碼、環境和 source／binary／log／產物身分。

## 剩餘範圍

Legacy 其餘矩陣／vector 回傳碼、完整配置與外部資產一致性、constructor 清理
仍需依 HPC-01C 補齊。本批不宣稱任意 PETSc 內部失敗都可恢復，也不保證
失敗 CLI 在退出前釋放每個舊有 PETSc resource。Output 仍需 root 集中場資料。
Sequential coupling 入口、其餘工具／CUDA／scheduler 與 HDF5 能力矩陣依
HPC-01D／06／09 追蹤；本批沒有 GPU、跨節點、效能或擴展性驗收。
