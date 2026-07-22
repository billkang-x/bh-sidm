# Baryon-aided Time-dependent SIDM Accretion Research Plan

日期：2026-07-04

## 项目目标

研究重子物理如何改变高红移黑洞种子在自相互作用暗物质（SIDM）晕中的时变吸积过程。目标是在现有 1D 球对称 SIDM 流体框架基础上，加入重子势场、重子气体吸积、冷却/反馈的最小模型，判断重子辅助是否能显著扩大轻种子黑洞成长为中等质量或超大质量黑洞种子的参数空间。

核心问题：

1. 重子势场是否能压缩 SIDM 内区并提高暗物质吸积率？
2. Eddington 限制的重子吸积是否能先把轻种子推入高效 dark Bondi/SIDM 吸积阶段？
3. 热传导、重子冷却、辐射反馈三者竞争时，黑洞增长是否仍能达到 Little Red Dots 所需的质量范围？
4. 结果对 halo mass、redshift、concentration、SIDM cross section、baryon fraction 和 feedback efficiency 有多敏感？

## 文献锚点

当前本地文献：

- `2607.02151v1.pdf`：Spherically Symmetric Fluid Simulations of Black Hole Accretion in Self-Interacting Dark Matter Halos。它提供了 1D 球对称 SIDM 流体、热传导、黑洞内边界吸积的基础框架。

关键外部文献：

- Feng, Yu, Zhong 2025, Dark Bondi Accretion Aided by Baryons and the Origin of JWST Little Red Dots: https://arxiv.org/abs/2506.17641
- Zhong, Yang, Yu 2023, The impact of baryonic potentials on the gravothermal evolution of self-interacting dark matter haloes: https://arxiv.org/abs/2306.08028
- Jiang et al. 2025, Formation of the Little Red Dots from the Core-collapse of Self-interacting Dark Matter Halos: https://arxiv.org/abs/2503.23710
- Feng, Yu, Zhong 2021, Seeding Supermassive Black Holes with Self-Interacting Dark Matter: https://arxiv.org/abs/2010.15132
- Outmezguine et al. 2022, Universal gravothermal evolution of isolated self-interacting dark matter halos for velocity-dependent cross sections: https://arxiv.org/abs/2204.06568
- Silverman et al. 2026, Mergers Matter: Gravothermal Collapse in Dwarf Halos with Self-Interacting Dark Matter: https://arxiv.org/abs/2606.02566

## 工作假设

先从最小可发表模型开始，而不是直接做完整多维辐射流体：

1. SIDM 仍使用 1D 球对称自引力流体加热传导模型。
2. 重子先作为外加且可随时间演化的球对称势场进入动量方程。
3. 黑洞质量增长分为两路：Eddington 限制的重子吸积和内边界 SIDM 通量吸积。
4. 重子冷却与反馈先用参数化模型，后续再升级为显式气体方程。
5. 优先研究高红移小晕到中等质量晕：`M_halo ~ 10^6-10^9 M_sun`，`z ~ 10-30`。

## 阶段计划

### 阶段 1：文献与理论框架

目标：把现有 SIDM-BH、baryon-aided dark Bondi、baryonic potential gravothermal evolution 三条线合并成一个统一模型。

交付物：

- 文献矩阵：模型假设、变量、方程、参数范围、主要结论。
- 模型设计文档：应该修改哪些方程、哪些量作为输入、哪些量作为输出。
- 最小可行科学问题：重子辅助是否改变 `M_BH(t)` 的数量级。

### 阶段 2：基线模型复现

目标：复现本地论文的无重子结果，至少复现趋势和量级。

交付物：

- 1D SIDM 流体求解器或等效简化版本。
- SIS 与 NFW baseline run。
- 热传导开/关对比。
- 边界、网格、时间步敏感性检查。

### 阶段 3：加入重子势场

目标：先研究静态或准静态重子势场如何改变 SIDM 密度、速度弥散和黑洞暗吸积率。

模型路线：

- 在动量方程中加入 `M_b(r,t)`。
- 测试 Hernquist、Plummer、isothermal gas core 或 compact exponential baryon profile。
- 参数扫描 `M_b/M_halo`、`r_b/r_s`、形成时间 `t_b`。

交付物：

- 有/无重子势场的 `M_BH(t)` 对比。
- 重子集中度对 dark accretion 的增强因子。
- `rho_DM`、`v_DM`、`u_DM` 的径向演化图。

### 阶段 4：加入重子吸积与反馈

目标：让黑洞先通过重子 Eddington accretion 成长，再触发更强 SIDM/dark Bondi accretion。

模型路线：

- `dM_BH/dt = dM_DM/dt + dM_b/dt`
- `dM_b/dt = min(dotM_Bondi,b, dotM_Edd)` 或参数化 Eddington ratio。
- 反馈先用中心加热/降低 baryon concentration 的参数化项。

交付物：

- 重子吸积、暗物质吸积分别贡献的质量历史。
- 判断何时进入 dark-accretion dominated 阶段。
- 与 Little Red Dots 所需 `M_BH ~ 10^6-10^8 M_sun` 的时间尺度比较。

### 阶段 5：参数空间与可观测连接

目标：把单个模拟升级为科学结论。

参数：

- `M_halo`、`z`、`c`、`M_seed`
- `sigma/m` 或 velocity-dependent `sigma(v)/m`
- `f_b`、`r_b`、baryon formation time
- Eddington ratio、duty cycle、feedback efficiency

交付物：

- 成功/失败相图：哪些参数能在给定红移前达到目标质量。
- 与 LRD black-hole mass、host mass、overmassive ratio 的比较。
- 对 SIDM cross section 的反推约束。

### 阶段 6：论文结构

初步论文题目：

`Baryon-aided Time-dependent SIDM Accretion onto High-redshift Black-hole Seeds`

论文骨架：

1. Introduction：JWST/LRD tension、SIDM seeding、为什么需要时变流体模型。
2. Model：SIDM conducting fluid + baryonic potential/accretion。
3. Numerical Method：operator splitting、inner boundary、validation。
4. Results I：baseline and baryonic potential。
5. Results II：baryonic accretion/feedback and parameter space。
6. Discussion：适用性、限制、观测预测、未来多维/宇宙学模拟。

## 风险与优先级

最高风险：

- 球对称可能高估吸积。
- 重子反馈可能破坏中心高密度环境。
- 内边界处理可能控制黑洞增长率。
- 常数大截面可能难以同时满足低红移约束。

优先解决：

1. 基线复现和内边界敏感性。
2. 重子势场是否带来数量级增强。
3. 速度依赖截面是否保留高红移快速增长。
4. 反馈是否把增强效应抹掉。

## 每轮进展汇报规则

后续每轮对话开始时，先汇报：

1. 当前阶段。
2. 已完成内容。
3. 最新发现或阻塞点。
4. 下一步将执行的具体任务。

