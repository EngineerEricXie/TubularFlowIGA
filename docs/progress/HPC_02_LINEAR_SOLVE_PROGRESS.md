# HPC-02 固定幾何 Newton 線性求解修正

本項支援 HPC-02 的原生數值驗收，不代表 HPC-04 的可配置求解器、
預條件器矩陣或網格擴展工作已完成。

## 原始失敗

`outputs/hpc02/volume/fixed-newton4/` 在原生 depth-4、27-cell、兩個
flow controllers 與 pressure gauge 案例完成七次 Newton 更新，但最後的
真實線性相對殘差 `2.5883782137720832e-4` 超過既有 `1e-10` gate，
因此退出 134。非線性殘差通過不能替代這個線性驗收。

## 隔離矩陣診斷

`outputs/hpc02/linear-probe/` 封存相同 fixture 第一次組裝的 PETSc
Jacobian 與 RHS。診斷直接包含既有測試中的幾何、Options 與初始場定義，
透過僅測試建置可用的 hook 複製矩陣；沒有更換模型或 RHS。
擷取及 replay 的原始 GMRES 數值完全一致，RHS norm 均為
`1.0945723152615692`。`comparison.json` 綁定來源、執行檔、矩陣與原始日誌。

環境為 PETSc 3.15.5、GCC 11.4、Open MPI 4.1.2、FP64、32-bit PetscInt。
各候選使用新的 KSP、相同 restart 30、LU shift `1e-20`、KSP rtol `1e-16`、
atol `1e-50`、divtol `1e5` 與最大 2,000 次迭代；left 候選使用
preconditioned norm，其餘為 unpreconditioned norm。所有真實殘差均另以
`MatMult`、`VecAXPY`、`VecNorm` 計算。

| 候選 | 迭代數 | 真實相對殘差 | 原生線性 gate |
|---|---:|---:|---|
| 原始右預條件 GMRES | 30 | 1.3822381914851506e-5 | 未通過 |
| GMRES，CGS refinement always | 5 | 6.9064944345501689e-6 | 未通過 |
| GMRES，modified Gram–Schmidt | 29 | 1.8073171434319128e-5 | 未通過 |
| FGMRES，右預條件 | 30 | 2.0549559498475232e-15 | 通過 |
| GMRES，左預條件 | 34 | 2.7502748312506434e-10 | 未通過 |

五項均回報 KSP reason 2，再次顯示遞迴殘差與正的收斂原因不能替代真實殘差。
這是單一擷取線性系統的證據，尚非完整非線性或 FSI 驗收。

## 修正與後續驗收

浸入式暫態 runtime 改用 FGMRES，保留既有右預條件 LU、restart、容許值及
迭代上限；固定 solver configuration hash 也改記 FGMRES。
FGMRES 保存預條件後的方向向量。根據本次相同矩陣的比較，推論此案例的
差異來自右預條件求解結果的建構方式；不能推廣為所有 GMRES 問題均需替換。
算法與方向保存方式可參考 [PETSc FGMRES 文件](https://petsc.org/release/manualpages/KSP/KSPFGMRES/)
及 [FGMRES 原始碼](https://petsc.org/release/src/ksp/ksp/impls/gmres/fgmres/fgmres.c.html)。
連結為目前文件，表格則來自本機 3.15.5 實測。

`fixed-newton-fgmres4/` 已完成完整、未放寬 gate 的 Newton 驗收，退出 0。
四個 threads 綁定 CPU 8／10／12／14，timeout 1,500 s。七次更新後
非線性殘差 `5.7228617089619486e-14`，最後真實線性相對殘差
`1.3514175383586189e-15`；全部七次線性求解的真實相對殘差均小於
`2.06e-15`。原生收斂、controller／gauge、iteration cap 與 accounting
gate 均通過，`accepted-evidence.json` 另行核對全部七次線性結果。

來源與執行檔另行封存；process wall 為 1,231.504684 s、peak RSS 為
106,364,928 bytes，14 次 assembly 為 1,218.242073 s、七次 setup 為
0.903129 s、linear solve 為 0.088257 s。這輪與其他正確性測試同時執行，
不能與先前 GMRES 失敗執行的時間比較來宣稱性能收益。

`fsi-fgmres4/` 已取得 accepted collection，原生 FSI 收斂、force／moment、
moving mass、wall leakage 及 continuity gates 均通過。與 HPC-00D 參考
的九場比較全部通過：速度相對 L2 `2.7784424817471747e-13`、壓力
`6.961917358512401e-12`，其餘場最大為 traction 的 `1.2213742720409472e-14`。
物理 fixture input identity 不變；修正後來源與執行檔已封存於 collection
目錄。該次 process wall 為 308.142277 s、peak RSS 為 73,314,304 bytes，
與 Newton 同時執行，不能作隔離效能證據。`fsi-fgmres2/` 與
`fsi-fgmres8/` 也已取得 accepted collection，合計 27 個場比較通過；
2／4／8 threads 的九個輸出場檔案雜湊逐一相同。每份 collection 都另行
封存來源、執行檔與基準比較。

修正後的 1-thread 完整場回歸也已通過。隔離矩陣
`outputs/hpc02/volume/isolated-matrix-fgmres/` 已完成八次原生執行與 72 個
場比較，最大 relative L2 為 `6.961917358512401e-12`。4-thread 端到端
速度比為 2.300008，RSS ratio 為 1.021657，FGMRES 額外向量成本已包含在
量測中；來源、binary 及原始證據已封存並重新核對，見
[FSI 效能驗收](HPC_02_FSI_PERFORMANCE.md)。無 OpenMP binary 也已重新
建置並通過九個完整場比較，並非以 OpenMP build 的單 thread 取代。
