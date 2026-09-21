# 人工正方體雙樹：多步被動氧氣示蹤

這是**人工、零耗氧的被動濃度示蹤功能案例**，不是肝臟生理氧合或組織供氧驗證。它延續[細薄壁人工雙樹 FSI 案例](CUBE_DUAL_TREE_FSI.md)的動脈、固定 Darcy 組織與靜脈五區網格；輸送期間將**已通過的準穩態 FSI 回饋後**血流／壁位置及 Darcy RT0 面流量固定，不在每個濃度步重解流體或壁。血管壁不透物質；交換僅發生於明示的四個動脈末端與四個靜脈末端。組織孔隙率人工設為 1、消耗率為 0；沒有血紅素結合、氧分壓、耗氧反應或病例材料資料。

專案自主 P1 tetra 輸送 FEM 以 backward Euler、單調圖擴散，在三區分別求濃度。動脈末端的**求解物質外流 `mol/s`**依原來 matching face 的水流量及相鄰組織 cell 體積，形成逐 cell 組織來源 `mol/(m³·s)`；組織使用已恢復且保守的 Darcy RT0 面通量對流。每個靜脈末端的組織物質外流除以該端 `m³/s` 水流量，得到對應的靜脈入口 donor 濃度。每步都檢查三區求解收斂原因、非負濃度、逐區帳本、兩次界面轉移與全域物質帳本；任一失敗不宣告案例完成。PETSc 只供稀疏代數與 MPI，不替代弱式／元素／組裝。

## 版本化輸入與 WSL 重跑

[`cases/idealized_cube_dual_tree_oxygen.json`](../cases/idealized_cube_dual_tree_oxygen.json) 明列人工入口 `1 mol/m³`、三區共用擴散係數 `10⁻⁶ m²/s`、`Δt=5 s`、最少 10 步、最多 60 步，以及「靜脈出口平均濃度達入口的 1%」停止 gate。這些只是讓傳播可見的數值測試值；達標時間**不是**生理通過時間。若尚未有細薄壁 flow run，先按[FSI 文件](CUBE_DUAL_TREE_FSI.md)產生。其後執行：

```bash
make -C solvers/cpu native_tet_cube_oxygen PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
python3 scripts/run_idealized_cube_dual_tree_oxygen.py /tmp/idealized-cube-refined-fsi-r8 /tmp/idealized-cube-oxygen --case cases/idealized_cube_dual_tree_oxygen.json --ranks 8
```

第二行會重新產生／分割幾何，核對其 manifest 與 flow run 的 hash，再從已有 PVTU 擷取帶全域 ID 的 FSI 後速度、位移及 Darcy RT0 面流量。若已持有 hash 相符的分區網格，可加 `--split-directory /tmp/idealized-cube-refined-submeshes` 略去重建；輸出仍會複製網格到自己的目錄。輸出目錄必須不存在，不會覆寫。WSL 執行 MPI／ParaView 需要允許本地 loopback socket。

## 在 ParaView 看組織網格與濃度前緣

開啟輸出目錄的 `transport/tissue_oxygen.pvd`，按 **Apply** 後選 **Coloring → concentration_mol_m3**，以 **Clip** 或 **Slice** 剖開 30 mm 正方體；Representation 選 **Surface With Edges** 就能看到**固定組織內部 tetra 網格**。按 **Play** 可看 0、5、10…秒的濃度變化，建議將 color range 固定，不要讓每幀自動重縮色階。另可開 `transport/artery_oxygen.pvd` 與 `transport/vein_oxygen.pvd`，用相同色階；三個 PVD 的時間軸相同。Darcy 流速本身仍在 flow run 的 `paraview/tissue_darcy_velocity.vtu`，可疊加 Clip／Glyph；這裡的濃度圖不是速度圖。

組織 mesh 座標在全部步數均**不變形**；動／靜脈輸送使用 FSI 後固定幾何。因此 Play 看到的是濃度前進，不是新的血管壁 FSI 時間動畫。檔案 `ledger.csv` 記錄每步三區庫存、兩端交換、靜脈出口濃度及收支缺陷；`provenance.json` 綁定案例、flow 摘要、網格、擷取場與求解器 SHA-256。

## 本地功能驗收

2026-09-20 WSL 的細薄壁 flow run 上，8-rank 與 2-rank 輸送各自完成 12 步（含 `t=0` 共 13 個狀態），在 `t=60 s` 的靜脈出口平均濃度約 `0.01071475877 mol/m³`，即入口的 `1.0715%`。最終較嚴格的守恆 gate 下，8-rank `ledger.csv` 的最大絕對全域收支缺陷約 `5.84×10⁻²⁰ mol/s`；ParaView 已讀入三區 PVD 的全部 13 個時間值與組織每步 101,716 tetra。組織最後一步節點濃度約 `0.00155–0.342 mol/m³`，在固定色階的剖面可看出空間變化。獨立 2-rank 從原人工幾何重新產網格，案例／flow／分區／擷取場／求解器 hash 與 8-rank 相同，39 個濃度場陣列按全域節點 ID 比對最大絕對差 `9.9067×10⁻¹² mol/m³`。另以最多 2 步、要求 90% 出口濃度的負向試驗，在尚未突破時明確 exit 2，沒有發布成功的 PVD 時間索引。可重跑比較：

```bash
python3 scripts/compare_idealized_cube_oxygen.py /tmp/idealized-cube-oxygen-verified-r2 /tmp/idealized-cube-oxygen-verified-r8
python3 -m unittest scripts.tests.test_idealized_cube_oxygen
make -C solvers/cpu native-tet-moving-species-transport-test
```

此結果證明專案自主輸送與明示來源／匯的**數值功能**，不證明真實組織可得到多少氧氣。入口濃度、孔隙率、擴散與零消耗的人工設定尤其不能推為肝臟參數。若未來要模擬可用氧，需另訂血氧形式／溶解與結合關係、壁交換、組織孔隙與耗氧模型、材料資料及相應實驗驗證；不會從此功能案例自行推定。
