# HPC-00D：浸入式與 FSI 完整數值參考

狀態：浸入式與 FSI 各兩次完整輸出、比較及拒絕測試已通過。
本批接續 [物種收支驗收](HPC_00D_TRANSPORT_BUDGET.md)，整項驗收見
[HPC-00D 報告](HPC_00D_REPORT.md)。

## 輸出與比較契約

兩個原生 fixture 增加可選 `--reference-output NEW_DIRECTORY`。完成原有
FD、Newton、守恆、牽引力、強耦合與狀態提交檢查後，才輸出完整數值場。
未指定此參數時沿用原本執行路徑；沒有修改生產數值核心。

[ReferenceStateOutput.hpp](../../solvers/cpu/tests/ReferenceStateOutput.hpp)
按物理量分檔，使用十進位整數 ID 與足以精確往返 double 的文字精度。
每份檔案記錄單位、列數、欄數及 SHA-256。最後才發布 `manifest.json`；
既有目錄會被拒絕，任何寫入錯誤都不發布完成 manifest。
這是小型數值參考輸出，不是 checkpoint 或大型分散式 I/O。

| 路徑 | 比較物理量 | ID |
|---|---|---|
| 浸入式 | 流體速度、壓力、流量控制器壓力 | 背景節點 ID、控制器 boundary label |
| FSI 流體 | 速度、壓力 | 背景節點 ID |
| FSI 膜與介面 | 膜位移、膜速度、牽引力、consistent nodal force | 不變的膜 global node ID |
| FSI 材料表面 | 位置、位移、速度 | 原始 source material vertex index |

所選案例都有 pressure outlet，因此沒有 gauge unknown；浸入式有一個
flow-rate controller，FSI 的兩個 port 都是 pressure control。
輸出器會輸出實際存在的 controller／gauge；若物理設定變更導致欄位集合
改變，比較器要求重新檢視案例與政策，不默默略過額外自由度。
浸入式輸出是已收斂的靜態 trial；FSI 是時間 `1.05 s`、step 1 的
accepted macro-step。兩者均不是完整續跑狀態。

[hpc_reference_states.json](../../benchmarks/hpc_reference_states.json) 在原生
參考執行前登記每個物理量的相對 L2 `1e-6` 與零參考絕對 L2 `1e-12`。
各物理量分別檢查，不合併壓力、速度、力與位移的範數；不調整 pressure
gauge，也不以大尺度座標代替微小材料位移的比較。所有原生物理 gate 保留。

[hpc_reference_states.py](../../scripts/hpc_reference_states.py) 的 `collect`
要求原生成功退出、完整 profile／原生診斷、有效場 manifest 與執行前後
不變的輸入／code／binary 雜湊。`compare` 核對兩個已接受 collection、
案例定義與輸入身分、場 schema／epoch／ID／單位／雜湊，再逐物理量比較。
只檢查 manifest 的 `native_gates_passed` 欄位不構成原生驗收證據。

## 命令與環境

```bash
make -C solvers/cpu reference_state_output_test
./solvers/cpu/reference_state_output_test outputs/hpc00/reference-states/writer-unit1
make -C solvers/cpu immersed_aneurysm_jacobian_test compliant_channel_fsi_test \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real -j2
python3 -m unittest discover -s scripts/tests -p 'test_hpc*.py' -v

python3 scripts/hpc_reference_states.py collect --case fsi \
  --output-dir outputs/hpc00/reference-states/fsi-reference --cpu 0 --timeout 1500
python3 scripts/hpc_reference_states.py collect --case fsi \
  --output-dir outputs/hpc00/reference-states/fsi-repeat --cpu 2 --timeout 1500
python3 scripts/hpc_reference_states.py collect --case immersed \
  --case-dir outputs/hpc00/serial-matrix/immersed-inputs \
  --output-dir outputs/hpc00/reference-states/immersed-reference --cpu 4 --timeout 300
python3 scripts/hpc_reference_states.py collect --case immersed \
  --case-dir outputs/hpc00/serial-matrix/immersed-inputs \
  --output-dir outputs/hpc00/reference-states/immersed-repeat --cpu 4 --timeout 300

python3 scripts/hpc_reference_states.py compare \
  outputs/hpc00/reference-states/immersed-reference outputs/hpc00/reference-states/immersed-repeat
python3 scripts/hpc_reference_states.py compare \
  outputs/hpc00/reference-states/fsi-reference outputs/hpc00/reference-states/fsi-repeat
```

重跑須使用新的 output 目錄。工作站為 WSL／i9-14900KF，GCC 11.4、
PETSc 3.15.5、OpenMPI 4.1.2，OpenMP／BLAS threads 設為 1。
`lscpu -p=CPU,CORE` 確認 logical CPU 0／2／4 分屬實體 core 0／1／2。
兩個 FSI 與浸入式驗證程序有重疊執行；此批時間不納入 HPC-00C 的
隔離效能矩陣。HPC-00C 先前已完成的數字仍屬先前記錄的執行檔版本。

本批確切來源與執行檔封存於
`outputs/hpc00/reference-states/reference-source-and-binaries.tar.gz`，共 138
個 source／build／binary 檔案。各 collection 的 `inputs.json` 保存每個
檔案的 SHA-256；封存不包含 PETSc 等外部工具鏈的安裝。
原始輸出與測試日誌亦位於 `outputs/hpc00/reference-states/`，不提交生成產物。

## 已完成驗收與剩餘工作

22 項原生 writer 檢查通過，涵蓋 `2^53` 以上 ID 的精確往返、數值精度、
完成標記時機、重複／負 ID、非有限值、I/O 失敗與既有目錄保護。
新增 11 項 Python 比較測試，包含微小位移、零參考、ID 不一致、截斷、
錯誤單位／epoch、缺失欄位與完成 manifest；HPC Python 測試合計 41 項通過。

浸入式兩次都通過完整九組 FD 與原生數值 gate，並輸出 216 個速度節點、
216 個壓力節點及一個控制器。三個物理量的相對 L2 全為 0。
另以實際輸出副本擾動壓力係數並更新測試副本的完整性 metadata，CLI 的
數值比較退出 `2`，壓力相對 L2 `2.02172704`；此副本明確標為故障注入，
不當作真實求解結果。紀錄在 `immersed-comparison.json` 及 `corrupt-pressure-test.json`。

FSI 兩次均正常退出，九個物理量的相對 L2 全為 0：速度／壓力各 216
個背景節點，膜位移／速度／牽引力／nodal force 各 9 個材料節點，
位置／位移／速度各 26 個原始表面頂點。四個強耦合迭代、原有 moving
mass／wall leakage／continuity 及合力／力矩 gate 均通過。
兩次子程序時間分別為 `630.261762 s`、`629.931064 s`，峰值 RSS 為
`75,780,096`、`72,523,776 bytes`；僅記錄本批資源使用，不據此宣稱加速。

七項原生拒絕測試亦通過：FSI／浸入式既有目錄、缺少參數、FSI 重複
參數、兩 rank 輸出，以及截斷 surface 的原生 gate 失敗，均退出 `1`。
失敗 fixture 沒有完成狀態 manifest，collector 保持 failed。
完整命令與退出碼在 `native-rejection-tests.json`。

最後重新核對四份 collection 的 profile／原生診斷／schema／source／binary、
兩份完整比較結果及封存內每個對應來源的 SHA-256。驗收彙整為
`outputs/hpc00/reference-states/accepted-evidence.json`，包含原始產物雜湊。
參考可用於下一階段數值比較；尚未執行 OpenMP、多 rank immersed／FSI
或跨節點平行化驗收。
