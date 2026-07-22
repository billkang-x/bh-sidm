# Research Progress Log

项目：Baryon-aided Time-dependent SIDM Accretion onto High-redshift Black-hole Seeds

## 2026-07-04

当前阶段：阶段 1，文献与理论框架。

已完成：

- 确认工作目录只有一篇本地文献：`2607.02151v1.pdf`。
- 阅读并归纳本地文献的核心内容：1D 球对称 SIDM 自引力流体、热传导、黑洞内边界吸积；SIS 环境中 `100 M_sun` 种子可在 2 Myr 内增长到约 `10^4 M_sun`，NFW 环境增长明显较弱。
- 确定研究方向：在该时变流体框架上加入重子辅助机制。
- 创建持续研究目标：后续每轮对话都先汇报研究进展。
- 保存研究计划到 `research_plan.md`。

最新判断：

- 最有价值的切入点是先加入重子势场，再加入 Eddington 限制的重子吸积，最后再考虑反馈。
- 这一路线能把本地论文的时变模拟优势与 2025 年 baryon-aided dark Bondi accretion 的科学问题连接起来。

下一步：

- 完成文献矩阵，明确每篇关键论文提供的模型组件。
- 写出最小可行模型的方程修改清单。
- 决定第一版计算是“静态重子势场”还是“时间增长的重子势场”。

补充进展：

- 已保存 `model_blueprint.md`。
- 决定第一版计算采用“静态 Hernquist 重子势场”，随后升级为时间增长的重子势场。
- 明确了最小方程修改：在 SIDM 动量方程中加入 `M_b(r,t)`，黑洞质量增长后续拆成 `dotM_DM + dotM_b`。
- 设计了前五组实验：baseline、static baryonic potential、growing baryonic potential、baryon+dark accretion、feedback robustness。
- 已创建 Python 原型结构：`sidm_bh/`、`tests/`、`configs/`。
- 已实现 `HernquistBaryons` 重子势场模块，包括 enclosed mass、density、potential、acceleration 和指数增长质量模型。
- 已加入 Eddington accretion rate 工具函数。
- 已保存第一组静态重子势场实验配置：`configs/experiment_1_static_baryons.yaml`。

下一轮起手任务：

- 运行并检查 `tests/test_baryons.py`。
- 根据测试结果修正重子模块。
- 开始整理 SIDM 基线求解器的数据结构和无量纲单位转换。

验证：

- `python -m unittest tests.test_baryons -v` 已通过，4 个测试全部 OK。

更新后的下一步：

- 整理 SIDM 基线求解器的数据结构。
- 从本地论文提取无量纲单位转换和初始条件生成方法。
- 将 `HernquistBaryons.acceleration_cgs()` 接入未来的动量方程源项。

## 2026-07-04 继续推进

当前阶段：阶段 2，基线模型复现的代码准备。

已完成：

- 新增 `sidm_bh/units.py`：实现论文式 (5)、(6) 的 `SimulationScales`，包括 `M0`、`v0`、`t0`、`p0`、`sigma/m` 标度和热传导标度。
- 新增 `sidm_bh/sidm.py`：实现 SMFP/LMFP 插值热传导、平均自由程、引力尺度高度和 Knudsen 数。
- 新增 `sidm_bh/halos.py`：实现 NFW 与 SIS profile 的密度和包围质量。
- 新增 `sidm_bh/mesh.py`：实现球对称有限体积网格，体积和界面面积与本地论文离散形式一致。
- 新增 `sidm_bh/state.py`：实现 primitive variables 与 conservative variables 的转换。
- 新增 `sidm_bh/initial_conditions.py`：实现 profile 采样和 hydrostatic pressure 初始条件生成。
- 修正 hydrostatic 外边界处理：`outer_pressure_code` 现在表示外界面压力，最外层 cell center 压强保持正定。
- 新增 baseline 复现实验配置：`configs/experiment_0_baseline.yaml`。

验证：

- `python -m unittest discover -v` 已通过，13 个测试全部 OK。
- 关键验证包括：物理单位热传导无量纲化后等于论文 Eq. (9)、NFW 质量构造正确、SIS 质量-密度关系正确、hydrostatic 初始条件可恢复等温 `rho ~ r^-2` 解。

最新判断：

- 代码已经具备搭建 baseline 1D solver 的基础积木。
- 下一步应实现最小 operator-splitting 框架中的源项与通量接口，然后先做无热传导/静态 halo 的 sanity check，再接 Roe flux 和 implicit conduction。

下一步：

- 写 `sources.py`：组合 `M_DM + M_BH + M_b` 的引力源项，并把 Hernquist 重子势场接入。
- 写 `accretion.py`：实现内边界 SIDM accretion rate 的 code/cgs 转换。
- 开始实现最小 hyperbolic flux 接口，为后续 Roe solver 做准备。

## 2026-07-04 继续推进 2

当前阶段：阶段 2，基线模型复现的数值框架搭建。

已完成：

- 新增 `sidm_bh/sources.py`：实现 `M_DM + M_BH + M_b` 的总包围质量、Hernquist 重子包围质量采样、无量纲引力加速度和球对称 Euler 源项。
- 源项形式与本地论文保守形式一致：`S = [0, -rho M/r^2 + 2p/r, -rho M u/r^2]`。
- 新增 `sidm_bh/accretion.py`：实现内边界 SIDM 吸积率 `dotM = -r_min^2 rho_min u_min`，并支持 code units 与 `M_sun/Myr` 互转。
- 新增 `sidm_bh/fluxes.py`：实现最小 convective Euler flux、`gamma = 5/3` 声速和最大信号速度接口。
- 扩展 `FluidState.from_conservative()`，支持 conservative variables 反变换回 primitive variables。
- 新增 `tests/test_sources_accretion_flux.py`，覆盖源项、重子质量采样、吸积率、单位转换、通量和 conservative roundtrip。

验证：

- `python -m unittest discover -v` 已通过，24 个测试全部 OK。

最新判断：

- 基线 solver 的基础接口已经比较完整：网格、状态、profile、单位、热传导 closure、源项、吸积率、通量都已经有测试覆盖。
- 还没有实现真正的时间推进；下一步需要搭建 CFL 时间步和球坐标 finite-volume update 骨架，然后才能接 Roe solver 和 implicit conduction。

下一步：

- 写 `time_integration.py` 或 `solver.py`：实现 CFL 时间步估计。
- 实现 spherical finite-volume convective update skeleton：`U^{n+1} = U^n + dt(S - div F)`。
- 先使用简单 Rusanov/Lax-Friedrichs flux 做 sanity check，再替换或扩展到 Roe flux。

## 2026-07-11 继续推进

当前阶段：阶段 2，基线模型复现的双曲时间推进已经可运行。

已完成：

- 扩展 `sidm_bh/fluxes.py`：新增 Rusanov（local Lax-Friedrichs）界面通量。左右态相同时严格退化为物理 Euler 通量。
- 新增 `sidm_bh/solver.py`：实现局部 CFL 时间步、零梯度边界、当前 SIDM 密度的包围质量积分，以及论文式 (14) 的球对称有限体积显式推进。
- 自引力质量积分使用无量纲关系 `dM/dr = r^2 rho`，并允许显式传入计算域内边界以内的残余暗物质质量。
- 双曲推进支持可选黑洞质量、重子包围质量和外部给定的暗物质包围质量，因此 baseline 与静态 Hernquist 势场可共用同一接口。
- 新增 `tests/test_solver.py`，覆盖 Rusanov 一致性、零梯度界面状态、局部 CFL、自引力质量积分、有限体积全域质量守恒，以及内边界流体损失与黑洞吸积诊断的一致性。
- 新增 `sidm_bh/conduction.py`：按论文式 (15) 和附录 B 式 (B2)-(B5) 实现 frozen-coefficient 隐式热传导、零梯度无热流边界和 Thomas 三对角求解。
- 新增组合稳定步长 `stable_timestep_code()`。除声波 CFL 外，它还限制显式引力在一步内相对于声速的速度改变量；这是初态突然打开中心黑洞引力时保持内能正定所必需的。
- 新增 `tests/test_conduction.py`，覆盖三对角求解、常温解不变、无通量边界热能守恒和温差衰减。

验证：

- `python -m unittest discover -v` 已通过，37 个测试全部 OK。
- `python -m compileall -q sidm_bh tests` 已通过。
- 使用论文 NFW 参数 `rs = 30 pc`、`rho_s = 3.7 M_sun/pc^3`、`M_BH = 100 M_sun` 和 256 个对数网格，对有/无热传导两条分支各进行 100 个 operator-split 步；演化后密度与速度弥散保持正定。
- 该短时测试中传导分支最内层 `v^2` 比绝热分支低约 4.1%，仅作为机制方向检查，不作为物理结果。
- 进一步进行 1000 个带传导组合步的稳定性测试，推进至约 `4.74e-3 Myr`，状态继续保持正定。

最新判断与限制：

- 当前已经闭合一阶 Rusanov 双曲推进与隐式热传导，但仍不是论文完整算法：尚未加入 MC 重构、Roe solver、黑洞质量逐步更新和长时输出驱动器。
- NFW 组合冒烟测试只推进约 `5.0e-4 Myr`，不能据此宣称复现论文 2 Myr 的质量增长曲线。
- 只使用声波 CFL 在强中心引力下会产生非正内能；引力源项步长约束必须保留，并在后续长时运行中监测是否过度限制效率。
- 在进行重子参数扫描前，必须先完成黑洞质量联动和 baseline 长时验证。

下一步：

- 增加无黑洞 NFW 热传导验证，优先对照论文 Figure 1 的 SMFP/LMFP 趋势。
- 增加完整演化驱动器，在每步用内边界质量通量同步更新 `M_BH`，并跟踪流体质量、黑洞质量和边界通量预算。
- 在 Rusanov 基线稳定后，再加入 MC 重构和 Roe flux，评估数值耗散对吸积率的影响。

## 2026-07-11 继续推进 2

当前阶段：阶段 2，完整 operator-splitting 演化与质量预算已经闭合，并开始进行论文 Figure 1 的 NFW 验证。

已完成：

- 新增 `sidm_bh/evolution.py`：组合稳定时间步、双曲更新、可选隐式传导和黑洞质量更新，形成可运行的演化驱动器。
- 每个接受步使用内边界 SIDM 通量同步更新 `M_BH`，并记录 `M_BH(t)`、`dotM_DM(t)`、网格内流体质量、外边界带符号通量、内边界回流和质量预算残差。
- 新增 `tests/test_evolution.py`：验证流体损失与黑洞增长闭合、外边界流出预算、内边界回流不扣减黑洞质量、可选传导子步和最大步数失败机制。
- 将单元中心电导率计算向量化，并增加与标量 Eq. (9) closure 的逐点一致性测试。

验证：

- `python -m unittest discover -v` 已通过，43 个测试全部 OK。
- `python -m compileall -q sidm_bh tests` 已通过。
- 论文 NFW+BH 参数下推进 `0.005 Myr`：绝热与传导分支均由 `100 M_sun` 增长到约 `100.299 M_sun`，峰值 SIDM 吸积率约 `79 M_sun/Myr`，相对质量预算残差小于 `9e-16`。
- Figure 1 无黑洞验证采用正确的 `r_min = 0.1 kpc`。在固定 `200 pc` 处，SMFP 的传导/绝热密度比在 128、256、512 网格分别为约 `0.717`、`0.718`、`0.725`，传导的额外中心稀释基本收敛。
- LMFP 首个快照的传导效应较弱，传导/绝热密度比随 128、256、512 网格从约 `0.944`、`0.954` 到 `0.964`。
- 128 网格推进到论文最晚快照：SMFP `47 Myr` 和 LMFP `754.92 Myr` 均显示中心密度显著降低及热结构重分布，质量预算残差约 `1e-15`，演化方向与 Figure 1 一致。
- 向量化后，128 网格 SMFP `47 Myr` 测试约为 `1.7e3 step/s`。

最新判断与限制：

- 黑洞质量与流体质量现在在离散层面严格闭合，长时预算可审计。
- Figure 1 的“传导相对绝热效应”呈收敛趋势，但相对于初态的绝对剖面仍有明显分辨率依赖，说明一阶 Rusanov 数值耗散不可忽略。
- 当前结果只能作为趋势验证；在复现论文定量曲线和开展重子增强因子扫描前，需要降低双曲通量的数值耗散。

下一步：

- 实现论文使用的 MC slope reconstruction，并先与 Rusanov 通量组合，分离“空间重构”与“Riemann solver”两部分误差。
- 随后实现 Roe flux，对 NFW Figure 1 和 NFW+BH 吸积率进行 128/256/512 网格收敛比较。
- 定量基线通过后，再运行静态 Hernquist 重子势场的第一组有/无重子对照。

## 2026-07-11 继续推进 3

当前阶段：阶段 2，MC 空间重构、Roe 通量和新的网格收敛验证已经完成。

已完成：

- 新增 `sidm_bh/reconstruction.py`：实现向量化 minmod、适用于非均匀对数网格的 MC slope，以及 `rho/u/p` primitive interface reconstruction。
- MC 使用真实 cell-center 距离计算左右与中心斜率，并使用 center-to-interface 距离外推；在非均匀网格的线性数据上保持精确，在非单调区域不产生新极值。
- 在 `sidm_bh/fluxes.py` 中实现论文附录 A7-A16 等价形式的 Roe flux，包括 Roe velocity/enthalpy、三条特征波与声学 Harten entropy fix。
- `solver.py` 和 `evolution.py` 现在支持显式选择 `constant|mc` reconstruction 与 `rusanov|roe` Riemann solver。
- 增加可关闭的 Roe-to-Rusanov 正定回退。所有本轮收敛实验均关闭回退，确认结果来自纯 MC-Roe。
- 新增 `tests/test_reconstruction_roe.py`，覆盖非均匀网格线性精度、单调性、静止接触间断、全超声速迎风极限和 MC-Roe 质量守恒。
- 新增 `convergence_results.md`，保存三种格式、三档网格以及 SMFP/LMFP 与 NFW+BH 的完整收敛表。

验证：

- `python -m unittest discover -v` 已通过，50 个测试全部 OK。
- `python -m compileall -q sidm_bh tests` 已通过。
- MC-Roe 在 NFW+BH 短时运行中关闭正定回退仍可稳定完成，质量预算残差保持在 `10^-15` 量级。
- Figure 1 SMFP 的 `rho_cond/rho_init` 在 MC-Roe 128/256/512 网格分别为约 `0.7230/0.7249/0.7257`，256 到 512 变化约 `0.11%`。
- Figure 1 LMFP 对应结果约为 `0.9274/0.9313/0.9329`，256 到 512 变化约 `0.17%`。
- NFW+BH 在 `0.02 Myr` 的 MC-Roe 黑洞质量为 `101.8708/101.8805/101.8828 M_sun`，256 到 512 变化约 `0.0023%`；峰值吸积率收敛到约 `119.8 M_sun/Myr`。

最新判断与限制：

- 此前绝对剖面的主要分辨率漂移来自 piecewise-constant reconstruction；MC-Rusanov 已消除大部分误差，Roe 在此基础上提供较小但可测的修正。
- 256 个对数网格已足以支持早期 NFW 定量运行，512 网格适合作为敏感性检查。
- 当前 MC 重构选择正定 primitive variables，而论文文字使用 state vector `U`；Roe 还加入了论文未明确写出的熵修正。两者是有意的稳健性选择，应在方法部分明确说明。
- 尚未完成论文 Figure 2/4 的完整 `2 Myr` SIS/NFW 质量曲线复现。

下一步：

- 使用 MC-Roe、256/512 网格运行论文 SIS 与 NFW 的 `2 Myr` baseline，比较 `M_BH(t)` 和峰值/晚期吸积率。
- 对 inner boundary、CFL 和 entropy-fix 参数做敏感性检查。
- baseline 定量通过后，进入静态 Hernquist 重子势场的首组增强因子实验。

## 2026-07-11 继续推进 4

当前阶段：阶段 2 的完整 2 Myr baseline 已完成，内边界、CFL 和 Roe 熵修正敏感性已量化。

已完成：

- 为 SIS 增加自然无量纲标度 `SimulationScales.for_singular_isothermal_sphere()`：取 `r0 = 1 pc`、`rho0 = rho_SIS(r0)`，因此 `v0 = sqrt(2)c_s`。
- 增加严格恒温 SIS 初态生成器，初始一维速度弥散在整个计算域保持 `c_s = 4.2 km/s`，不受有限外边界零压近似影响。
- 增加 `sidm_bh/fast_evolution.py`：Numba 编译的 MC-Roe、显式源项、隐式传导和黑洞质量联动长时内核，可在指定物理时刻保存完整径向状态。
- 增加 fast/reference 交叉验证。短时 NFW 运行中两者步数和黑洞质量一致，密度最大相对差小于 `8e-16`、速度弥散差约 `1.5e-14`。
- 增加精确冻结引力 kick 作为可选源项积分；短时结果与严格 Euler 一致，但长时步数最终仍由双曲 CFL 控制，因此正式 baseline 继续使用论文的 Euler 源项形式。
- 完成 NFW 256/512 网格和 SIS 2 Myr 基线，并完成 inner boundary、CFL、entropy fix 敏感性矩阵。
- 新增 `baseline_2myr_results.md`，保存完整质量历史、吸积率与敏感性表。

验证：

- 全量测试现为 53 个，fast/reference 一致性测试通过。
- NFW 256 网格：`M_BH(2 Myr) = 129.0756 M_sun`，峰值/末期率约 `124.27/2.84 M_sun/Myr`；512 网格为 `129.1406 M_sun`，与论文 Figure 4 的约 `130 M_sun` 一致。
- SIS nominal `r_min = 0.001 pc`、128 网格：`M_BH(2 Myr) = 9646.8 M_sun`，与论文 Figure 2 的约 `10^4 M_sun` 一致；增长主要发生在前 0.8 Myr。
- NFW 内边界减半/加倍得到 `123.79/134.36 M_sun`，相对 baseline 约变化 `-4.1%/+4.1%`。
- NFW CFL `0.1-0.4` 与 entropy fix `0-0.2` 对最终质量影响低于 `0.001%`。
- SIS `r_min = 0.001/0.002/0.005/0.01 pc` 得到 `9646.8/10304.9/11699.1/13292.3 M_sun`；累计质量保持 `10^4 M_sun` 量级，但归一化明显受吸收半径控制。
- SIS 在固定 `r_min = 0.005 pc` 下，128 到 256 网格变化约 `1.1%`，CFL `0.2-0.8` 变化约 `0.012%`，entropy fix `0-0.2` 变化约 `0.02%`。
- 所有完成运行的相对质量预算残差保持在约 `1e-15` 到 `1e-13`。

最新判断与限制：

- MC-Roe 已定量复现论文 SIS/NFW 的 2 Myr 质量端点。
- CFL 与 Roe 熵修正不是主要不确定性；inner absorbing boundary 是主导系统误差，SIS 尤其明显。
- SIS nominal 固定 `0.001 pc` 边界需要约 607 万步，即使 JIT 后仍显著昂贵；峰值吸积率随边界变化很大，不应作为稳健物理预测。
- 后续重子增强必须使用完全匹配的网格和内边界做有/无重子比值，并至少报告两个 `r_min`，否则可能把边界效应误判为重子增强。

下一步：

- 进入阶段 3，运行静态 Hernquist 重子势场的首组匹配对照。
- 首先使用计算可控的 NFW 256 网格，并在 `r_min = 0.0025/0.005/0.01 pc` 下比较 `M_BH(t)` 和增强因子。
- SIS 重子实验先使用 `r_min >= 0.002 pc`，nominal 极小边界仅做少量确认运行。

## 2026-07-11 阶段 2 收尾

当前阶段：阶段 2 已完成；下一轮进入阶段 3。

本轮完成：

- 新增 `scripts/run_stage2_case.py`，可逐案例运行并立即保存带物理单位的质量历史、吸积率历史和完整径向快照。
- 新增 `scripts/plot_stage2_results.py`，生成 heat-flow on/off 增长曲线及六时刻 `rho/u/v` 径向剖面图。
- 完成 NFW 与 SIS nominal `2 Myr` heat-on/off 对照，共保存 201 个时刻。
- 完成外边界减半/加倍检查，并按计算域对数跨度调整网格数以保持近似相同的 `Delta ln r`。
- 新增 `stage2_validation_report.md`，逐项核对阶段 2 验收标准并给出完成决定。

核心结果：

- NFW heat-on/off 的 `M_BH(2 Myr)` 为 `129.08/234.36 M_sun`，末期率为 `2.84/57.00 M_sun/Myr`。
- SIS heat-on/off 的 `M_BH(2 Myr)` 为 `9646.8/32528.3 M_sun`，末期率为 `155.7/15886 M_sun/Myr`。
- 关闭传导使最终质量提高 `1.82` 倍（NFW）和 `3.37` 倍（SIS），最终内区密度提高约 `25.7` 倍和 `104` 倍。
- NFW `r_max = 2500/5000/10000 pc` 的最终质量跨度约 `1.7e-6`；SIS 控制设置 `r_max = 500/1000/2000 pc` 的跨度约 `1.4e-4`。
- 因此外边界影响可忽略，内吸收边界仍是主导数值系统误差。

交付文件：

- `results/stage2/*.npz`：基线与外边界原始数据。
- `results/stage2/figures/*.png`：增长与径向剖面图。
- `stage2_validation_report.md`：阶段 2 正式验证报告。
- `baseline_2myr_results.md`、`convergence_results.md`：详细敏感性与收敛数据。

阶段判断：

- 1D solver、SIS/NFW baseline、heat-flow on/off、网格/内外边界/CFL/熵修正敏感性和可复查输出均已完成。
- 阶段 2 正式结束。
- 下一步进入阶段 3：静态 Hernquist 重子势场的匹配有/无重子对照。

## 2026-07-11 阶段 3 启动

当前阶段：阶段 3，静态 Hernquist 重子势场；首批匹配矩阵已完成。

已完成：

- 按 `paracloud-skill` 核对 HPC 安全流程。当前主矩阵本地 Numba 可完成，未请求或保存任何凭据。
- 新增 `sidm_bh/stage3.py`，统一构造 `DM + BH + static Hernquist` 势下的匹配静水平衡初态。
- 新增 `scripts/run_stage3_static_case.py`，逐案例保存 2 Myr 质量历史、吸积率、径向状态、固定重子包围质量与元数据。
- 增加固定重子势下 fast/reference 一致性测试和初态测试。
- 完成 `f_b = 0.01/0.05/0.16` 与 `a_b/r_s = 0.001/0.01/0.1` 的 3x3 主矩阵，以及无重子匹配控制。
- 对一个增强案例和一个抑制案例完成 `r_min = 0.0025/0.005/0.01 pc` 三边界验证。
- 新增 `scripts/analyze_stage3_static.py`、CSV 汇总、三张结果图和 `stage3_static_baryon_report.md`。

主要发现：

- 匹配无重子控制在 2 Myr 内吸积 `22.12 M_sun` 暗物质。
- 紧致重子势 `a_b/r_s = 0.001` 的暗增长增强为约 `29.8/61.3/99.2` 倍，对应递增的三个重子分数。
- 中等集中度 `a_b/r_s = 0.01` 的增强为约 `12.3/30.0/41.8` 倍。
- 扩展重子势 `a_b/r_s = 0.1` 不增强吸积，而是抑制为控制的约 `0.88/0.70/0.50`。
- 代表增强案例的边界因子为 `22.3/30.0/43.5`，代表抑制案例为 `0.77/0.70/0.65`；边界改变幅度但不改变符号。

最新判断与限制：

- 阶段 3 的首个关键结果不是“重子总是增强”，而是增强存在明显的集中度阈值。
- 静态平衡协议排除了瞬时打开势场的冲击，但尚未描述重子质量的形成历史。
- `a_b/r_s = 0.01` 的高重子分数案例在 2 Myr 时吸积率仍上升，需要更长时间或增长势模型解释。
- 完整 3x3x3 边界矩阵再叠加截面变化适合迁移到 Paracloud；届时需用户通过会话环境提供凭据，绝不落盘。

下一步：

- 比较静态平衡势与逐渐组装的 Hernquist 势，判断增强是否依赖初态协议。
- 在紧致度阈值附近加密 `a_b/r_s`，定位增强到抑制的转变。
- 准备 Paracloud Slurm 参数矩阵脚本；获得凭据后运行完整边界和截面扫描。

## 2026-07-11 阶段 3：重子势组装历史

当前阶段：阶段 3；静态 Hernquist 终态势与有限时间组装势的代表性对照已完成。

已完成：

- 在参考求解器和 Numba MC-Roe 求解器中加入有限时间 smoothstep 重子质量组装；默认静态路径保持不变。
- 新增时变外势 fast/reference 守恒变量一致性测试，以及 `scripts/run_stage3_assembly_case.py` 单算例驱动器。
- 对 `f_b = 0.05`、`a_b/r_s = 0.01/0.1` 完成 `T_asm = 0/0.05/0.2/0.5/1.0 Myr` 的 2 Myr 扫描。
- 对 `T_asm = 0.5 Myr` 的两个代表点完成 `r_min = 0.0025/0.005/0.01 pc` 三边界复核。
- 新增 `scripts/analyze_stage3_assembly.py`、两份 CSV、四张结果图和 `stage3_assembly_history_report.md`。

主要发现：

- 紧致案例的响应对组装时间非单调；`T_asm = 0.5 Myr` 在 2 Myr 内吸积 `1258.23 M_sun` 暗物质，是静态平衡协议的 `1.90` 倍。
- 扩展案例发生定性翻转：静态平衡为无重子控制的 `0.70`，组装势则为 `2.16-3.23`，说明此前的“抑制”依赖预先热压支撑的初态协议。
- 在扩展案例 `0.1 pc` 处，静态平衡与 `T_asm = 0.5 Myr` 的终态密度分别为无重子的 `6.30` 和 `15.89` 倍，径向结构支持组装诱发收缩的解释。
- 扩展案例的组装增强在三个内边界为 `2.45/2.59/2.80`，符号翻转不受所测内边界改变。

最新判断与限制：

- 重子集中度不是唯一控制量；重子形成时间与 SIDM 热/动力学时间的相对大小同样重要。
- 静态平衡矩阵描述终态势下的平衡响应，不能替代形成历史实验。
- 当前 smoothstep 是受控模型而非宇宙学重子增长史；紧致案例在 2 Myr 时仍未饱和。

下一步：

- 在 `a_b/r_s` 与 `T_asm` 二维平面加密，定位非单调峰值和静态/组装差异的转变区。
- 增加 `sigma/m` 变化，比较重子组装时间与 SIDM 传导时间的耦合。
- 将完整浓度-组装时间-内边界-截面矩阵迁移到 Paracloud Slurm。

## 2026-07-12 阶段 3：HPC 三维参数面

当前阶段：阶段 3；重子集中度、组装时间和 SIDM 截面的三维扫描已完成。

已完成：

- 使用独立远程工作目录部署项目；通过离线轮子安装 Numba/llvmlite，未在任何文件中保存登录凭据。
- 生成 164 项严格清单：4 个截面匹配无重子控制，以及 `8` 个集中度、`5` 个组装时间、`4` 个截面的 160 个有重子算例。
- 主 Slurm 作业 `40734527` 按 QOS 限制打包为 8 个 worker，全部 `COMPLETED/0:0`；边界作业 `40734548` 的 5 个 worker 同样全部完成。
- 下载并验证全部 164 个主结果和 10 个边界结果；HPC 与本地 `sigma/m = 50` 控制质量严格一致。
- 新增 `scripts/analyze_stage3_hpc_matrix.py`、`scripts/analyze_stage3_hpc_boundary.py`、三份 CSV、JSON 摘要、四张图和 `stage3_hpc_parameter_scan_report.md`。

主要发现：

- 160 个组装势算例全部高于截面匹配无重子控制，最弱案例仍为 `1.138` 倍；静态平衡协议中的抑制未在组装参数面恢复。
- 全局最大为 `sigma/m = 10 cm^2/g`、`a_b/r_s = 0.01`、`T_asm = 0.5 Myr`，2 Myr 内吸积 `3422.10 M_sun` 暗物质，增强 `102.28` 倍。
- 紧致势存在有限组装时间最优值；`sigma/m = 10-50` 时约为 `0.5 Myr`，`sigma/m = 100` 时在最紧致点移到 `1 Myr`。
- `a_b/r_s >= 0.07` 时瞬时组装最强，延迟组装随 `T_asm` 增加而减弱。
- 三个关键点的边界因子分别为 `99.76/102.28/103.38`、`20.61/26.20/36.43` 和 `1.1388/1.1380/1.1413`。

最新判断与限制：

- 集中度是主控制量，但紧致势下组装时间与 SIDM 截面发生明显耦合。
- 最优 `T_asm` 随截面移动与时间尺度匹配假说一致，但尚未直接计算传导/动力学时间，不能视为已证明机制。
- 当前三维矩阵固定 `f_b = 0.05`、单一 NFW 晕和 2 Myr 终点。

下一步：

- 计算供给半径附近的局部动力学时间、碰撞/传导时间和组装时间比值，检验非单调峰值机制。
- 机制确认后再扩展重子分数与暗晕质量，避免盲目增加参数维度。

## 2026-07-12 阶段 3：时间尺度与传导机制

当前阶段：阶段 3；紧致重子势的时间尺度诊断、组装时间加密和 heat-off 判别已经完成。

已完成：

- 新增 `sidm_bh/timescales.py`，按演化器同一质量与电导率闭合计算动力学、碰撞、传导、流入时间和 Knudsen 数。
- 新增按向内质量通量定义的对数半径中位供给半径，以及 3 项基础单元测试。
- 对 160 个有重子 HPC 算例在组装中点、组装终点和 2 Myr 终点生成 Hernquist 半径与实际供给半径诊断。
- 回到 `t = 0` 无重子初态，在未来供给半径重新计算时间尺度，避免用已经演化的热状态循环解释最优时间。
- 将紧致势组装时间扩展到 `2 Myr`，并在峰值附近加密；补充 heat-off 完整组装时间曲线。
- 新增 `scripts/analyze_stage3_timescales.py`、`scripts/analyze_stage3_mechanism_followup.py`、多份 CSV/JSON、四张机制图和 `stage3_timescale_mechanism_report.md`。

主要发现：

- 四个有传导紧致最优点的 `T_opt/t_dyn(r_feed,t=0)` 为 `0.985-1.163`，对数 RMS 仅 `0.046 dex`。
- `T_opt/t_cond` 的 RMS 为 `0.131 dex`，仍为同阶但分散更大；`T_opt/t_coll` 跨 `0.739-7.422`，不能预测最优时间。
- 加密后的离散最优组装时间为 `0.60/0.50/0.65/0.75 Myr`，对应 `sigma/m = 10/30/50/100 cm^2/g`。
- heat-off 分支在 `1.0 Myr` 达到约 `2.94` 倍增强，说明动力学本身产生宽最优窗口，但无法产生有传导分支的 `27-102` 倍放大。
- 所有传导最优点在 Hernquist 半径和供给半径均为向外热流；传导移走压缩热并维持收缩，而不是向中心加热。

最新判断与限制：

- 当前证据支持“动力学时间确定供给范围，传导决定放大幅度”的组合机制，而不是碰撞时间共振。
- `sigma/m = 10` 的 `0.5-0.7 Myr` 结果仅相差约千分之一，应称为宽匹配平台而非尖锐共振。
- 供给半径仍从演化结果提取，尚未形成完全独立的事前预测公式。

下一步：

- 在 `f_b = 0.01/0.16` 和不同暗晕质量下，用统一的早期供给半径定义检验 `T_opt ~ t_dyn` 缩放。
- 若缩放保持，再将其整理为阶段 3 的预测模型并正式收尾。

## 2026-07-12 阶段 3：重子分数与暗晕质量缩放

当前阶段：阶段 3；事前预测缩放矩阵及峰值加密已经完成。

本轮完成：

- 定义固定特征密度、固定浓度的自相似 NFW 暗晕族，以
  `M_halo=1e6 M_sun`、`rho_s=3.7 M_sun/pc^3`、`r_s=30 pc` 为锚点，
  对应浓度 `c=3.9201`。
- 将结果反推的供给半径替换为事前定义 `r_supply=3.5 a_b`，并在无重子
  初态计算 `T_pred=t_dyn(r_supply,t=0)`。
- 在 `f_b=0.01/0.05/0.16`、`M_halo=1e6/1e7/1e8 M_sun`、
  `a_b/r_s=0.01`、`sigma/m=50` 下完成 39 个主矩阵和 21 个峰值加密案例。
- 外边界随 `r_s` 缩放，网格从 256 增至 270/284 格以保持对数分辨率。
- 新增可复现清单、Slurm 脚本、分析表、四联图和
  `stage3_baryon_halo_scaling_report.md`。

主要发现：

- `T_pred` 在三种质量下为 `0.564-0.595 Myr`，符合固定密度自相似族预期。
- 未修正的动力学预测对全部 `f_b=0.05` 案例精确到 25%，但在低、高重子
  分数下发生方向一致的系统偏移。
- `f_b=0.01/0.05/0.16` 的峰值倍率中位数为 `0.380/1.197/1.984`，得到
  经验关系 `T_opt ~= 1.06 t_dyn(3.5 a_b) (f_b/0.05)^0.603`。
- 最优吸积暗物质质量对重子分数呈次线性缩放，指数为 `0.44-0.62`。
- 固定种子和固定物理内边界协议下，最优吸积量随暗晕质量缓慢下降；这是
  自相似性被固定物理尺度打破的效应，尚不能解释为宇宙学质量关系。

验证：

- Slurm 作业 `40734715` 和 `40734740` 全部以 `0:0` 完成，错误日志为空。
- 60 个结果文件齐全，最大质量预算残差为 `1.82e-12` code units。
- 本地完整测试为 63 项，全部通过。

下一步：

- 运行固定 `M_BH/M_halo` 与 `r_min/r_s` 的小型自相似伴随矩阵，分离固定
  种子/边界对暗晕质量趋势的贡献，然后再正式收尾阶段 3。

## 2026-07-12 阶段 3 正式收尾：自相似闭合检查

完成内容：

- 新增两组伴随协议：缩放 `M_BH/M_halo` 与 `r_min/r_s` 但固定物理截面；
  以及进一步缩放物理截面、保持无量纲输运完全相同的严格自相似协议。
- 完成 12 个有/无重子 HPC 案例，作业 `40734760` 全部 `0:0` 完成，错误
  日志为空。
- 新增 `similarity_summary.csv`、`similarity_statistics.json` 和
  `stage3_similarity.png`。

关键结论：

- 固定物理种子/边界/截面的有重子增长-质量斜率为 `-0.255`。
- 缩放种子和内边界后斜率变为 `+0.744`，证明原负斜率主要来自固定物理
  种子和吸收边界，而不是暗晕质量本身。
- 同时保持无量纲截面不变后，有重子与控制斜率均严格恢复为 `1.000`，
  三种质量的 `Delta M/M_halo` 相对展宽仅 `2.46e-14` 和 `1.53e-14`。
- 固定物理截面仍使增强因子随质量从 `57.68` 降至 `13.02`，说明质量族中
  剩余的非自相似部分来自 SIDM 输运强度的变化。

阶段判断：

- 用户指定的静态/准静态 Hernquist 重子势阶段 3 已正式完成。
- 结论不自动推广到 Plummer 或气体核 profile。
- 下一步正式进入阶段 4：动态重子吸积与反馈。

## 2026-07-12 阶段 4：动态重子吸积与势场膨胀反馈

当前阶段：阶段 4；第一轮无反馈基线与参数化反馈矩阵已完成。

实现内容：

- 在 Numba MC-Roe 演化器中加入有限 Hernquist 重子储库和 Eddington 上限
  吸积，并严格分离 `dM_DM/dt`、重子气体流入率与黑洞实际保留的重子质量。
- 采用 `dM_b,BH=(1-epsilon_r)dM_gas`，同时记录辐射质量等价、剩余储库、
  两路累计质量和两路瞬时吸积率。
- 加入参数化反馈
  `a_b(t)=a_b,0[1+E_fb/E_bind]^eta`，其中
  `E_fb=epsilon_f epsilon_r Delta M_gas c^2`，
  `E_bind=G M_b^2/(6a_b,0)`。
- 新增 `sidm_bh/stage4.py`、阶段 4 单案例驱动器、两组矩阵及峰值加密脚本、
  配置文件、分析表、统计 JSON、两张四联图和正式报告。

无反馈发现：

- `f_Edd=0` 严格回归阶段 3：最终质量差 `8.54e-7 M_sun`，最终密度最大
  相对差 `4.87e-10`。
- 满 Eddington 时，10/100/1000 `M_sun` 种子的暗物质增长占比分别为
  `98.99%/98.50%/96.92%`；在组装完成后的测量区间始终由暗物质主导。
- 每进入黑洞 1 `M_sun` 重子会额外催化 `1.76/1.06/0.37 M_sun` 暗物质。
- 反向作用同样显著：SIDM 增长把重子 Eddington 质量相对孤立种子放大
  `15.23/4.86/1.92` 倍，形成双向正反馈。
- 2 Myr 内没有案例达到 `1e4 M_sun`；最大终值为 `3545.13 M_sun`。

反馈发现：

- 重子催化反转阈值经加密后为 `epsilon_f ~= 2.04e-6` (`eta=0.5`) 和
  `1.08e-6` (`eta=1`)。
- 在首个明确反转点 `epsilon_f=3e-6`，反馈能量只有约 `0.098 E_bind`，
  势场尺度仅增大 `4.8%/9.7%`；紧致势增强对小幅膨胀非常敏感。
- 极端 `epsilon_f=0.01, eta=1` 将尺度半径增大 85 倍，暗物质和总增长仅
  保留无反馈值的 `2.49%/2.83%`。

限制与判断：

- 当前 Eddington 比是外部给定的，尚无 Bondi 供给、显式重子气体或冷却。
- 反馈是瞬时整体膨胀，且结合能只含 Hernquist 自束缚；反转阈值是模型
  参数，不能直接解释为真实反馈效率。
- 下一步应实现 `min(dotM_Bondi,b, f_Edd dotM_Edd)`，并把暗晕和黑洞对
  重子的束缚加入有效 `E_bind`，再开展 10-100 Myr 长时增长。

## 2026-07-13 阶段 4：Bondi-Eddington 转换与有效结合能

实现内容：

- 将重子供给升级为
  `dotM_gas=min(dotM_Bondi, f_Edd dotM_Edd)`，逐时记录 Bondi 上限、
  Eddington 上限、实际重子率和限制通道。
- 加入固定 `rho_b,inf`、`c_s,b`、`v_rel` 和 Bondi `alpha` 参数。
- 将反馈归一化结合能扩展为 Hernquist 自束缚、NFW 暗晕势和初始黑洞势
  三部分；基准总结合能为 `2.926e50 erg`，是原自束缚值的 `2.463` 倍。
- 完成 9 个气体密度-声速案例和 18 个 Bondi+有效结合能反馈案例。

主要发现：

- 9 个气体案例严格分成 3 个初始 Eddington 受限、3 个全程 Bondi 受限和
  3 个 SIDM 驱动的 Bondi-to-Eddington 转换。
- 三个转换发生在 `0.89-1.61 Myr`；模拟转换质量与解析临界质量只差
  `0.8%-1.3%`，确认转换由 SIDM 增长黑洞质量触发。
- 保留的重子质量跨 `0.173-19.800 M_sun`，而暗吸积只跨
  `1276.99-1297.78 M_sun`；气体状态主要控制直接重子通道。
- 有效结合能将反馈反转阈值提高到：转换气体轨道
  `4.54e-6/2.49e-6`，Eddington 饱和轨道 `4.77e-6/2.67e-6`
  （分别对应 `eta=0.5/1`）。
- 即使 `epsilon_f=3e-5`，反馈也只把转换时刻从 `0.89` 延迟到
  `0.90 Myr`，没有阻止 SIDM 将重子通道推入 Eddington 区域。

数值审计：

- 作业 `40738000` 的 9 个 Bondi 案例和 `40738013` 的 18 个有效反馈案例
  全部 `0:0` 完成，错误日志为空。
- 最大 SIDM 质量预算残差为 `5.87e-13`，本地测试增至 75 项且全部通过。
- 初次作业 `40737992` 因远端 NumPy 不支持新版梯形积分 API，在积分前
  统一失败且未生成结果；改为兼容实现后重跑成功，不影响科学数据。

下一步：

- 将 `rho_b,inf` 与 `c_s,b` 连接到演化中的 Hernquist 储库和反馈状态，
  使势场膨胀同时稀释/加热 Bondi 环境，再判断通道转换是否仍然存在。

## 2026-07-13 阶段 4：Hernquist-Bondi 环境闭环

实现内容：

- 令 Bondi 密度随剩余储库质量和 Hernquist 尺度半径按
  `rho_inf=rho_0(M_rem/M_0)(a_0/a)^3` 演化。
- 将反馈能量无重复地分成加热比例 `chi_heat` 和膨胀比例
  `1-chi_heat`，并按
  `c_s^2=c_s0^2+gamma(gamma-1)E_heat/M_rem` 更新声速。
- 在 MC-Roe 每个时间步重新计算动态 Bondi 上限，并输出密度、声速、
  限制通道及双交叉诊断。
- 完成 65 个 2 Myr 案例，包括冻结/动态环境匹配对照、三种能量分配、
  `epsilon_f=1e-7` 到 `1e-3` 扫描及阈值加密。

主要发现：

- 无反馈动态环境仍在 `0.89 Myr` 进入 Eddington；单纯储库消耗不足以
  改变原结论。
- 纯膨胀阻断阈值位于 `1.3e-4` 到 `1.7e-4`；首个阻断案例把气体密度
  降到初值的 `0.232`，并把暗吸积压低 `21.64%`。
- 50/50 混合阈值位于 `5e-5` 到 `7e-5`；纯加热阈值位于
  `4e-5` 到 `5e-5`。纯加热首个阻断点的声速为 `17.04 km/s`，暗吸积
  只降低 `0.093%`，说明气体限制器与暗势响应可以明显分离。
- 纯膨胀 `epsilon_f=1.3e-4` 在 `1.26 Myr` 进入 Eddington，随后于
  `1.97 Myr` 返回 Bondi，出现持续约 `0.71 Myr` 的暂态 Eddington 窗口。
- 动态 Bondi 会削弱自己的反馈燃料；在 `epsilon_f=1e-3` 强反馈端，
  动态环境比冻结供给保留更多暗吸积，因为气体先被稀释/加热，后续反馈
  能量随之下降。

数值审计与限制：

- 作业 `40738045/40738053/40738059/40738084` 给出全部 65 个有效结果；
  最大质量预算残差 `7.18e-13`，本地 79 项测试全部通过。
- 当前加热是无冷却的一区瞬时混合闭包，纯加热阈值可能偏低；下一步应
  加入有限冷却时间的泄漏热库，检查阻断阈值是否保持。

## 2026-07-13 阶段 4：有限冷却热库

实现内容：

- 将累积热能改为
  `dE_th/dt=q_heat dotM_gas-E_th/t_cool`，每个 MC-Roe 时间步用解析指数式
  更新，并用 `expm1` 保持长冷却时间极限的精度。
- 输出累计热注入、当前热库能量和冷却损失；省略冷却时间或显式使用无穷
  均严格回归无冷却模型。
- 对纯加热和 50/50 混合反馈各取临界与强反馈两档，完成 48 个
  `t_cool=0.001-100 Myr`、无冷却及阈值加密案例。

主要发现：

- 临界混合反馈 `epsilon_f=7e-5` 的阻断冷却时间为 `2-3 Myr`；强混合
  `1.4e-4` 降至 `0.15-0.20 Myr`。
- 临界纯加热 `epsilon_f=5e-5` 需要 `10-30 Myr`，说明上一轮最低纯加热
  阈值基本属于近绝热结果；强纯加热 `1e-4` 只需 `0.30-0.50 Myr`。
- 首个阻断点的热能保留比例从强混合的 `0.237` 到临界纯加热的
  `0.984`，反馈强度越接近无冷却阈值，对热泄漏越敏感。
- 多个案例出现短暂 Eddington 窗口后回到 Bondi；仅看 2 Myr 终态比值会
  错判是否真正阻断转换。
- 纯加热扫描中暗吸积几乎不变；混合反馈中延长热保持会更早切断 Bondi
  燃料，减少后续膨胀能注入，因此暗吸积略微回升。

数值审计与限制：

- 作业 `40738103/40738125/40738132/40738143` 给出全部 48 个有效结果，
  最大质量预算残差 `7.43e-13`，80 项本地测试全部通过。
- 初次无冷却清单暴露 CRLF 后缀导致的 `inf*0` NaN；异常文件未用于结论，
  已通过解析器正规化、`expm1` 核和独立重跑三重修正。
- 当前仍是给定 `t_cool` 的一区模型；下一步应使用随密度、温度和金属丰度
  变化的物理冷却函数，并评估混合反馈的势场再收缩。
# 2026-07-13 Stage 4: physical Cloudy cooling and trapping sensitivity

- Replaced prescribed cooling times with a density-, temperature-, and
  metallicity-dependent no-UV-background Grackle/Cloudy equilibrium table.
- Added a self-consistent temperature/MMW solve and an implicit nonlinear
  thermal-reservoir update to the MC-Roe evolution kernel.
- Completed 39 two-Myr cases spanning four feedback configurations, five
  metallicities, and cooling-rate suppression down to `1e-6`.
- Standard optically thin cooling defeats both pure-heating configurations;
  even `epsilon_f=1e-3` enters Eddington at `0.89-0.92 Myr`.
- Mixed feedback with `chi_heat=0.5` and `epsilon_f=1e-3` remains Bondi
  limited at every metallicity because `a_b` expands by about `1.9x` and the
  gas density falls by `6.8-7.0x`, despite retaining less than 1% of the
  injected thermal energy.
- Pure heating prevents the transition only when the cooling-rate multiplier
  is reduced from `1e-5` to `3e-6`, requiring roughly `1e5-3e5` suppression
  relative to the optically thin rate.
- Metallicity changes the thermal state modestly but never changes the
  limiter class in this matrix.
- Jobs `40738208` and `40738246` completed with `0:0` exits; 39 outputs are
  complete, the maximum mass residual is `7.34e-13`, and all 86 tests pass.
- The next discriminating test is reversible Hernquist-radius evolution with
  recontraction plus an optical-depth estimate, rather than a wider
  metallicity scan.
# 2026-07-15 Stage 5 entry: cosmological and observable feasibility

- Added a Planck matter-plus-Lambda cosmic-time budget and adopted the
  observation-aware LRD target bracket `1e5,1e6,1e7 M_sun`.
- Added `M200c + z + concentration` NFW construction, virial velocity, and a
  black-hole influence-radius resolution diagnostic.
- Found that the legacy fixed-density NFW anchor corresponds to `z=25.44` at
  its original concentration.
- Completed paired 20-case fixed-physical and fixed-dimensionless inner
  boundary matrices. Their mass trends reverse sign and differ by up to
  `621.9x`; unresolved fixed-seed halo-mass trends are not scientific.
- Completed a 20-case resolved matrix with `M_seed/M200=1e-4` and
  `r_min/r_influence=0.4167`. Dark matter supplies `92.8%-98.6%` of 2 Myr
  growth; five cases reach `1e5 M_sun`.
- A conservative 2 Myr dark burst followed only by baryonic Eddington growth
  converts 14 of 228 redshift-target combinations from infeasible to
  feasible.
- In a resolved `M200=1e9 M_sun`, `z=30`, `M_seed=1e5 M_sun` no-feedback
  run, the black hole reaches `1e7 M_sun` at `10.74 Myr` and
  `1.252e7 M_sun` at 12 Myr; dark matter supplies `91.8%` of added mass.
- Completed a 144-case transport screen. Every `M200-z` pair is optimized by
  `c=8` and `sigma/m=10 cm2/g`; the centers are in the SMFP regime, where
  increasing the cross section lowers conductivity and suppresses inflow.
- Completed an 81-case baryon-frontier screen. The compact, high-baryon
  corner (`f_b=0.16`, `a_b/r_s=0.003`) is optimal at all three redshifts;
  the best `z=30` case reaches `6.957e6 M_sun` in 2 Myr.
- A ten-case timing refinement gives an internal optimum at
  `T_asm=1.25 Myr` for both `z=20` and `z=30`. The latter reaches
  `7.008e6 M_sun`, with dark matter supplying `98.47%` of added mass; no
  refined 2 Myr case reaches `1e7 M_sun`.
- Completed a 34-case physical-cooling feedback matrix on the optimal
  frontier. At `epsilon_f=1e-3`, mixed feedback retains `99.825%-99.928%`
  of the `7.008e6 M_sun` ceiling; even `epsilon_f=1e-2` retains
  `99.249%-99.332%`.
- All low-density transition cases still enter Eddington by `0.58-0.69 Myr`.
  The compact baryon reservoir has `1.056e57 erg` effective binding energy,
  so the strongest mixed case expands it by only about `4.8%`.
- Pure-heating feedback is erased by Cloudy cooling and retains `99.997%` of
  the mass; pure expansion is more effective but still retains `98.712%` at
  `epsilon_f=1e-2`.
- Jobs `40739143/40739246/40739322/40758549/40758751/40758762/40759419/40759476/40761285`
  completed successfully with empty error logs; the latest matrix has a
  maximum mass residual of `2.42e-12`.
- Next: test the compact `a_b/r_s=0.003` frontier against a resolved feeding
  radius and reversible Hernquist-radius evolution before extending it in
  time.

# 2026-07-19 Stage 5: frontier closure and trusted threshold crossing

- Completed 26 one-axis closure cases, 72 interaction cases, 40 accepted
  general convergence cases, and 15 targeted compactness/convergence cases.
- The raw `c=12`, `sigma/m=100 cm2/g`, `a_b/r_s=0.0005` corner reaches a
  nominal `8.889e7 M_sun`, but fails inner-boundary convergence by 56.5% and
  is rejected.  The nominal `2.196e7 M_sun` resolved-scale high-sigma branch
  fails by 63.6%; moderate `sigma/m=3` branches fail by 12.4%-13.5%.
- The surviving `c=8`, `sigma/m=1 cm2/g` branch has an internal compactness
  maximum at `a_b/r_s=0.020`, rather than at the most compact reservoir.  Its
  nominal 2 Myr mass is `1.6850e7 M_sun`, with 98.81% of growth supplied by
  dark matter.
- Grid refinement changes the peak mass by only 0.298%.  The two smallest
  inner-boundary solutions differ by 7.14%, so the exact terminal mass does
  not yet satisfy the pre-declared 5% criterion.
- All peak variants cross `1e7 M_sun` within 1.52-1.76 Myr; the conservative
  completed boundary case ends at `1.4523e7 M_sun`.  The threshold crossing
  is robust even though the terminal mass remains boundary sensitive.
- Extending the current run merely to reach `1e7 M_sun` has no additional
  scientific meaning.  The next useful extension is a matched 3-5 Myr
  post-crossing saturation test after improving the feeding/capture boundary
  and reversible Hernquist response; a direct 10 Myr terminal mass would be
  model extrapolation rather than a robust prediction.
- Jobs `40761689/40761906/40762091/40764056/40764460/40765136/40765277`
  completed with `0:0` exits.  Job `40762942` yielded 40 valid cases; two
  optional ultra-small-boundary tasks reached the configured maximum steps.

# 2026-07-19 Stage 5: light-seed and inner-boundary closure

- Added a mass-conserving dark-matter feeding reservoir that separates fluid
  supply across `r_feed` from black-hole capture. Outside the influence
  region it uses the gamma=`5/3` Bondi coefficient; inside the calibrated
  region it switches to direct flux capture and drains the reservoir on a
  local free-fall time.
- The accepted influence-gated closure reproduces three resolved `1e5 M_sun`
  seed runs to better than `4e-6` in final mass. A pure local-Bondi closure
  was explicitly rejected because it is invalid once `r_feed` lies inside
  the influence radius.
- In the fixed `M200=1e9 M_sun`, `z=30` trusted halo, `1e2`, `1e3`, and
  `1e4 M_sun` seeds end at `103.3`, `1103`, and
  `1.33e4-1.72e4 M_sun`. The stellar-remnant cases pass 5% boundary and grid
  convergence and none reaches an LRD threshold.
- Light-seed failure is a capture bottleneck, not a supply failure. At the
  smallest feeding radius the three seeds receive `1.42e4`, `5.28e4`, and
  `6.74e6 M_sun` of dark supply but capture at most `0.20%`.
- At nominal `lambda=0.25`, `2e4 M_sun` fails at both boundaries,
  `3e4-3.75e4 M_sun` is boundary ambiguous, and `4e4 M_sun` is the lowest
  tested seed that crosses `1e7 M_sun` at both boundaries. The exact
  intermediate-seed threshold is not universal: `4e4 M_sun` fails for
  `lambda=0.20` and succeeds for `lambda=0.30`.
- The robust project claim is now rapid SIDM amplification of an existing
  intermediate/heavy seed. The optimized two-Myr mechanism does not solve
  the stellar-remnant light-seed problem.
- Jobs `40775904/40776966/40777024/40778274/40778815/40779138` completed all
  59 cases with `0:0` exits and empty error logs.

# 2026-07-20 Stage 5: velocity-dependent transport and applicability map

- Added the Born/Rutherford differential cross section with separate
  Maxwell-averaged viscosity factors `K3` and `K5`. The LMFP conductivity
  uses `sigma0 K3`, while the SMFP conductivity uses `(sigma0 K5)^-1`, rather
  than one virial-speed effective cross section.
- A logarithmic quadrature and 4097-point lookup table resolve the low-speed
  relative-velocity tail. The constant limit agrees with the original fast
  solver to `2e-10`; all 104 local tests and 12 targeted remote tests pass.
- Job `40782036` completed all 46 cross-section calibration tasks with `0:0`
  exits, empty error logs, and maximum mass residual `1.92e-11`.
- At the smaller calibration boundary, constant `sigma/m=1` reaches
  `1.560e7 M_sun`. Of the three pre-declared `sigma0/m=30` velocity brackets,
  `w=10` and `30 km/s` end at `1.816e6` and `4.901e6 M_sun`, while
  `w=100 km/s` reaches `1.065e7 M_sun` at 1.86 Myr.
- A virial-scale effective cross section near unity is not predictive by
  itself. The failed `w=30 km/s` model has virial `sigma0 K3/m=0.700` but
  initial inner `sigma0 K3/m=0.0248`; the successful `w=100 km/s` model has
  initial inner `sigma0 K3/m=0.875` and evolves strongly with time.
- All calibration terminal masses remain more than 5% boundary sensitive.
  The threshold classification survives for the constant control and the
  high-transport velocity model, but exact terminal masses are provisional.
- Built a 368-case six-parameter controlled map: four microphysics models,
  each with 28 one-axis points and 64 Latin-hypercube points over
  `M200=1e8-1e10 M_sun`, `z=15-30`, `c=4-12`, seeds `1e2-1e5 M_sun`,
  `a_b/r_s=0.005-0.08`, and `T_asm=0.25-1.75 Myr`.
- Job `40782909` is the first 46-point map chunk. Coordinator `40783158`
  automatically chains the remaining chunks, analysis, and an 80-case
  smaller-boundary/double-grid refinement selected around the target surface.
- The chained screen completed 366/368 points. Two high-mass/high-resolution
  points (`task_280` and `task_341`) reached the original `1.6e8` step cap,
  without numerical crashes; their recovery job `40795953` uses an
  `8e8`-step cap. Main analysis is held behind result monitor `40795957` and
  will resume automatically once both NPZ files are present.
- Recovery `40795953` completed both points with `0:0`; the final map now has
  `368/368` screen results and `80/80` boundary/grid refinement results.
  Final statistics are in `results/stage5/applicability_map_statistics.json`
  and `results/stage5/applicability_refinement_statistics.json`.
- The audited classifications are: constant `sigma/m=1` has 6 robust
  successes, 2 boundary-ambiguous points, and 8 failures; the high-transport
  Rutherford model (`sigma0/m=30`, `w=100 km/s`) has 2 robust successes, 7
  ambiguous points, and 7 failures; the low-transport model has 16 failures;
  the matched model has 15 failures and 1 ambiguous point.
- No audited point passes the full 5% terminal-mass boundary-plus-grid
  criterion. The final scientific classifier remains crossing of `1e7 M_sun`,
  not an exact 2 Myr terminal mass.
