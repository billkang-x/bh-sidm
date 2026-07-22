# Literature Notes

## Local baseline paper

### Meng et al. 2026

本地文件：`2607.02151v1.pdf`

主题：Spherically Symmetric Fluid Simulations of Black Hole Accretion in Self-Interacting Dark Matter Halos

用途：

- 提供 1D 球对称 SIDM conducting-fluid 方程。
- 使用自引力、热传导、operator splitting、Roe solver、implicit conduction。
- 黑洞增长由内边界质量通量给出。
- 结果显示 SIS 比 NFW 更容易产生快速黑洞增长。

可继承组件：

- SIDM hydrodynamic equations。
- SMFP/LMFP conductivity interpolation。
- Inner-boundary accretion prescription。
- Heat-flow on/off diagnostic。
- SIS/NFW baseline comparison。

需要扩展：

- 加入重子质量 `M_b(r,t)` 到引力项。
- 加入重子吸积 `dotM_b`。
- 加入反馈或冷却的参数化。
- 测试 velocity-dependent SIDM cross section。

## Baryon-aided dark accretion

### Feng, Yu, Zhong 2025

链接：https://arxiv.org/abs/2506.17641

主题：Dark Bondi Accretion Aided by Baryons and the Origin of JWST Little Red Dots

关键点：

- 提出重子先通过 Eddington accretion 增长黑洞，然后黑洞通过 dark Bondi accretion 继续快速增长。
- 目标是解释 `z ~ 4-11` Little Red Dots 中 `~10^7 M_sun` 黑洞。
- 该工作偏半解析，我们的机会是做时变 SIDM 流体版本。

可借鉴：

- baryon-aided growth scenario。
- LRD 目标质量与时间尺度。
- 暗物质吸积与重子吸积贡献分解。

## Baryonic potential in SIDM halos

### Zhong, Yang, Yu 2023

链接：https://arxiv.org/abs/2306.08028

主题：The impact of baryonic potentials on the gravothermal evolution of self-interacting dark matter haloes

关键点：

- 中心重子势场会显著改变 SIDM 晕的 gravothermal evolution。
- 重子势场可加速 expansion 和 collapse 阶段。
- 提供了把重子势场加入 SIDM gravothermal fluid model 的理论动机。

可借鉴：

- 静态 baryonic potential benchmark。
- baryon concentration 参数化方式。
- 与 N-body 校准的演化趋势。

## SIDM LRD population

### Jiang et al. 2025

链接：https://arxiv.org/abs/2503.23710

主题：Formation of the Little Red Dots from the Core-collapse of Self-interacting Dark Matter Halos

关键点：

- 用 Monte Carlo merger trees 研究 SIDM core collapse 产生 LRD 的统计可能性。
- 说明高浓度 halo 和合适 SIDM cross section 对 LRD abundance 很关键。

可借鉴：

- 后续连接 halo population 的参数范围。
- LRD mass function 作为观测检验。

## Early SIDM SMBH seeding

### Feng, Yu, Zhong 2021

链接：https://arxiv.org/abs/2010.15132

主题：Seeding Supermassive Black Holes with Self-Interacting Dark Matter

关键点：

- 重子可加速 SIDM gravothermal collapse。
- 高红移 SMBH 可能来自 SIDM halo 中心 collapse。

可借鉴：

- 高红移 halo 条件。
- baryon-accelerated collapse 的理论解释。

## Velocity-dependent SIDM evolution

### Outmezguine et al. 2022

链接：https://arxiv.org/abs/2204.06568

主题：Universal gravothermal evolution of isolated self-interacting dark matter halos for velocity-dependent cross sections

关键点：

- 研究 velocity-dependent cross section 下的 SIDM gravothermal evolution。
- 可为后续从常数截面升级到 `sigma(v)/m` 提供方法。

可借鉴：

- velocity-dependent scattering timescale。
- collapse-time rescaling。

## Assembly history risk

### Silverman et al. 2026

链接：https://arxiv.org/abs/2606.02566

主题：Mergers Matter: Gravothermal Collapse in Dwarf Halos with Self-Interacting Dark Matter

关键点：

- merger history 会改变 SIDM heat transport 和 core collapse。
- 安静并合历史更容易 collapse，持续并合可能阻止 collapse。

对本项目的意义：

- 1D isolated halo 模型需要在 Discussion 中承认 assembly history 风险。
- 后续 population-level 工作应加入 merger history。

