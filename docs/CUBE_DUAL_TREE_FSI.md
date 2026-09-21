# 人工正方體雙樹功能案例

此案例**不是肝臟解剖或生理驗證**。輸入是 [`cases/idealized_cube_dual_tree.json`](../cases/idealized_cube_dual_tree.json) 的 SI centerline 節點、逐段內半徑、固定壁厚、人工材料與邊界條件。兩棵樹各有七段、四個終端；動脈從正方體 `x=0` 面流入，靜脈從 `x=0.03 m` 面流出。各段以圓柱、分叉節點以球接合成封閉樹表面。Gmsh OCC Boolean 將正方體分成互斥的動脈腔、動脈壁、靜脈腔、靜脈壁與扣除兩棵樹外體積的固定 Darcy 組織，再 fragment 成共節點 tetra。Gmsh **只處理幾何與網格**；所有 FEM 弱式／元素／組裝在本專案實作，PETSc 僅做稀疏代數與 MPI。

本地需要 `python3` 的 Gmsh、NumPy；求解需要 MPI/PETSc。PNG 後處理另外需要 VTK Python bindings 與 Matplotlib；沒有它們仍可讀 `summary.json` 和以 ParaView 開啟 PVTU。

以這條**已求解的 FSI 後流場**做多步人工氧氣示蹤、查看組織 tetra 與靜脈出口濃度，見[人工雙樹氧氣功能案例](CUBE_DUAL_TREE_OXYGEN.md)。輸送採固定流場，不是隨氧濃度再耦合的暫態 FSI。

## WSL 重跑

```bash
make -C solvers/cpu native_tet_cube_dual_tree_fixed_flow PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
python3 -m unittest scripts.tests.test_idealized_cube_dual_tree
python3 scripts/run_idealized_cube_dual_tree_fixed_flow.py /tmp/idealized-cube-a --ranks 1
python3 scripts/run_idealized_cube_dual_tree_fixed_flow.py /tmp/idealized-cube-b --ranks 2
python3 scripts/compare_idealized_cube_dual_tree.py /tmp/idealized-cube-a /tmp/idealized-cube-b
python3 scripts/visualize_idealized_cube_dual_tree.py /tmp/idealized-cube-a
python3 scripts/export_idealized_cube_paraview.py /tmp/idealized-cube-a
```

較細薄壁版本（另用不存在的輸出路徑）：

```bash
python3 scripts/run_idealized_cube_dual_tree_fixed_flow.py /tmp/idealized-cube-refined --case cases/idealized_cube_dual_tree_refined.json --ranks 8
python3 scripts/visualize_idealized_cube_dual_tree.py /tmp/idealized-cube-refined --case cases/idealized_cube_dual_tree_refined.json --magnification 5000
python3 scripts/export_idealized_cube_paraview.py /tmp/idealized-cube-refined --magnification 5000
```

輸出路徑必須不存在；wrapper 會拒絕覆寫。WSL／容器需允許本地 MPI loopback socket。`summary.json` 記錄案例、求解器、幾何及分區 manifest hash、正向流量帳本、收斂原因、FSI 殘差及限制。`fields/` 有八組 rank-owned PVTU／VTU：初始動脈／Darcy／靜脈，兩壁位移，以及移動後的動脈／Darcy／靜脈。`overview.png` 是由 PVTU 重建的示意；位移為使形狀可見而放大，倍率與真實最大位移均標在圖上。

2026-09-20 WSL 最終程式的兩次獨立 1／2-rank 執行在 `/tmp/idealized-cube-fsi-final-r1`、`/tmp/idealized-cube-fsi-final-r2` 通過：70,791 tetra，最小 scaled Jacobian `0.0127652`；回饋後動脈→Darcy `6.544721406613861×10⁻⁹ m³/s`、Darcy→靜脈 `6.544721406613780×10⁻⁹ m³/s`、靜脈出口 `6.544721406613783×10⁻⁹ m³/s`。兩壁最大位移向量模分別 `2.8535×10⁻⁸`、`2.2939×10⁻⁸ m`；FSI 位移相對殘差 `2.3430×10⁻⁵`，最小流體變形 Jacobian `0.9999969`。兩次幾何／分區 hash 一致，36 個 PVTU 場陣列按全域 ID 比對最大絕對差 `3.2685×10⁻¹²`。數值只是人工參數的功能證據，尤其極小壁位移**不**證明生理合理。

### 較細網格與較薄壁的替代輸入

[`cases/idealized_cube_dual_tree_refined.json`](../cases/idealized_cube_dual_tree_refined.json) 保留相同的人工雙樹中心線、腔內半徑及功能測試材料／邊界條件；把壁厚由 `0.60 mm` 改為 `0.10 mm`，並以曲率網格選項把管口圓周約 `8` 個節點增加到 `32` 個。這只是較合理的**視覺／網格比例試值**，不是測得的血管壁厚度。最小 mesh size 設為 `0.04 mm`，使薄壁有足夠單元；細網格不等於分岔球—圓柱接合處已達切向連續。舊案例與其已通過證據保留，兩個輸入不得混稱同一次驗證。使用 `--case cases/idealized_cube_dual_tree_refined.json` 重跑上述 wrapper；繪圖也須指定相同的 `--case`。

本機 Gmsh 4.8.4 對 `0.10 mm` 動脈外樹沿舊反向建構順序時，Boolean 留下非零 containment 殘餘並產六個而非五個體積；此輸入顯式改用正向外樹建構，殘餘為零且五區互斥。驗收仍要求區域／界面及正 Jacobian 檢查，不以調大 Boolean 容差掩蓋重疊。

此輸入在本地幾何階段生成 `271,870` 個 tetra，最小 scaled Jacobian `0.00330819 > 0.001`；五個分區與 12 個 matching interface 標籤、原生 reader／P2 拓樸及正 Jacobian 都通過。這些僅是**網格驗收**，不能代替細網格上的血流／Darcy／FSI 收斂與守恆驗收。

本地 2／8-rank 各自從此輸入獨立重建並完成準穩態求解（約 `1504`／`575 s`，兩次有部分同時執行，**不是**嚴格的平行加速基準）：回饋後動脈→Darcy `8.133254228498714×10⁻⁹ m³/s`、Darcy→靜脈 `8.133254228498833×10⁻⁹ m³/s`、靜脈出口 `8.133254228498835×10⁻⁹ m³/s`；最大 Darcy cell 收支缺陷 `8.17×10⁻²³ m³/s`，FSI 位移殘差 `2.80×10⁻⁵`，最小流體變形 Jacobian `0.999979`，兩壁最大位移各 `6.92×10⁻⁹`／`5.15×10⁻⁹ m`。兩次案例／求解器／幾何／分區 hash 相同，36 個 rank-owned PVTU 場陣列按全域 ID 比較最大絕對差 `6.89×10⁻¹¹`。此值比厚壁舊網格的位移**更小**；改壁厚同時改變有限元素離散與牽引，不能把差異直接解釋為壁材料的生理趨勢。`paraview/` 已從 8-rank 結果匯出完整固定 Darcy 網格／RT0 速度、體積帳本及兩樹 `0/1` 狀態 PVD；使用顯式 `--magnification 5000`，一階 tetra 顯示 Jacobian 比皆為正（四區最小約 `0.67`），而 `×20000` 會倒置薄壁。ParaView `pvpython` 已開啟動／靜脈 pipe 與壁的四個 PVD，驗證各自兩個狀態均可載入。這仍僅是功能驗證，且未達 C2 壓力耦合及 C3 完整 FSI gate。

## 求解與驗收界線

- 本案例的人工入口速率為 `0.002 m/s`，血液測試密度／黏度為 `1000 kg/m³`／`0.004 Pa·s`；動脈四個交換端帽與靜脈出口設人工 `0 Pa`，Darcy 靜脈四端帽亦設 `0 Pa`，其餘組織外面不通量。Darcy mobility 為 `10⁻⁸ m²/(Pa·s)`。壁厚 `0.0006 m`、人工 Young modulus `2×10⁴ Pa`、Poisson ratio `0.3`；不能將這些值當成肝臟材料或病人 BC。
- 動／靜脈各用本專案 P2 速度／P1 壓力穩態 tetra 流體 FEM；固定組織用 P1 Darcy 與保守 RT0 面通量。四個動脈末端的 FEM 面通量依共面、唯一 owner 與法向映到相鄰 Darcy cell；四個靜脈末端的 Darcy 面流量作為靜脈入口流率控制。逐末端及全域 `m³/s` 平衡、正 Jacobian、PETSc 收斂原因與局部 Darcy cell deficit 均為拒絕 gate。
- 兩個血管壁各用本專案小應變 P1 tetra 線彈性；外壁和根部端帽固定，內壁受流體 Cauchy 牽引。固體位移透過彈性 ALE 延拓到流體網格，再重求流體與 Darcy，並以新流體牽引重算壁位移。程式要求位移相對固定點殘差不大於 `10⁻³`、移動後 Jacobian 為正與流量再次閉合。這是**準穩態、小應變的雙向流體—血管壁幾何／牽引功能測試**，不是時間準確的暫態 FSI。
- 目前水力跨域是**單向的保守流量傳遞**：動脈末端面通量送進 Darcy，Darcy 靜脈端面通量送進靜脈。尚未施加血管—組織壓力連續或雙向壓力／流量迭代；不能宣稱完整壓力耦合灌流。也尚未驗收 FSI 界面功率、GCL、暫態儲存及 trial／commit／rollback。組織網格始終固定。人工材料數字只用於測試，不對應肝臟或病人。

因此 C0／C1 幾何路線已單獨驗收；C4 的獨立重跑／跨 rank 全場比對子集也通過，但 C2／C3 尚缺完整 gate，不能將 C2–C4 勾成首版全部完成。其他既有心臟／瓣膜／心肌程式未刪除，但不屬於此案例。

## 圖與 ParaView 表面序列怎麼看

`overview.png` 左邊紅／藍線是**輸入**的兩棵 centerline，不是求解出的速度流線；綠色短箭頭是從固定組織求得的 RT0 Darcy 通量中取樣後**只顯示方向**，箭長不表示速度大小。右邊是血管壁 FEM 節點的 `x–y` 投影：灰色為變形前，紅／藍色為位移放大後的位置；第三維被投影掉，倍率與未放大的位移數值標在圖上。上方的 `m³/s` 數字是實際求解流量，不能把圖上箭頭密度或線條粗細當作流量大小。

ParaView 開啟 `paraview/arterial_pipe_x1000.pvd` 與 `paraview/venous_pipe_x1000.pvd`（或所指定的安全倍率），按 **Apply**、將兩者設為 **Surface**，再按動畫工具列的 **Play** 或切到 `0`／`1`；兩棵 pipe 的表面會隨狀態改變。`*_actual.pvd` 是**真實**求解位移，約 `10⁻⁸ m`（舊案例），以 30 mm 正方體全域視角幾乎看不出變化；`*_x1000.pvd` 僅把幾何位移放大方便觀察，`displacement_m` 陣列仍保留未放大的 SI 物理值。先前產出的 `×30000` 動脈腔有倒置 tetra，**不能當作可信的幾何變形圖**；匯出器現在以一階 tetra 的正 Jacobian 拒絕這類過度放大，但此 gate 不是完整高階曲面自相交證明。另有 `arterial_wall_*.pvd`、`venous_wall_*.pvd` 可看固體血管壁而非流體腔表面。按 **Coloring → displacement_m → Magnitude** 可查看位移分布；可用 **Reset Camera** 再縮放到血管附近。

PVD 的 `0`／`1` 表示「參考」／「準穩態 FSI 回饋後」兩個**狀態**，不是兩個物理秒數或完整暫態 FSI 動畫；目前沒有可宣稱的連續時間運動。匯出時會以雙精度重新建立高階 tetra 的表面點；先前 `×30000` PVD 雖可由 `pvpython` 載入，但其中動脈腔倒置，故載入成功**不代表幾何有效**。舊案例在 `×1000` 的一階 tetra 最小顯示 Jacobian 比約 `0.9969`，真實最大動脈位移為 `2.85350722669×10⁻⁸ m`。

### Darcy 內部速度與體積變化

同一 ParaView 輸出目錄還有 `tissue_darcy_velocity.vtu`：這是**固定組織的完整 54,978 個 tetra**，cell data `darcy_velocity_m_s` 是專案保守 RT0 重建的 cell-centroid 速度（`m/s`），`darcy_speed_m_s` 是其模長。組織是厚實的立方體，直接顯示外表面會遮住內部；在 ParaView 對此 `.vtu` 加 **Clip** 或 **Slice**，Representation 選 **Surface With Edges** 顯示內部 tetra，再以 `darcy_speed_m_s` 著色。可疊加 `tissue_darcy_direction_arrows.vtp`，其中約 180 個箭頭在實際 cell centroid，方向來自 RT0 速度，但箭長固定為 1 mm 以便顯示；不要從箭長推速度。實際速度大小看 cell array，或使用 **Cell Centers → Glyph** 自行調整箭頭密度。`pvpython` 已確認兩檔可載入。

`arterial_lumen_deformation.vtu`、`venous_lumen_deformation.vtu` 與兩個 `*_wall_deformation.vtu` 含每個 tetra 的 `deformation_jacobian` 及 `relative_cell_volume_change = J−1`；以 cell data 著色可看到非零體積應變。`volume_audit.json` 對一階 tetra 幾何的參考／變形後體積求和：動脈腔 `+0.00249095%`，靜脈腔 `+0.00193601%`；動脈壁 `−0.00138038%`，靜脈壁 `−0.00107052%`。每棵樹的腔體增加量與壁體減少量約相反，符合外壁固定的契約；組織網格座標前後完全一致，**正方體整體與 Darcy 組織不變形**。這些微小數值在原尺寸視角看不見，不代表只有位移場、沒有幾何變化；也不能把 30,000 倍顯示幾何拿來報物理體積變化。原生 PVTU 的實際座標已移動，對它再直接套一次 Warp By Vector 會重複位移。
