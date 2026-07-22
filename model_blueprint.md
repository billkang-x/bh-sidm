# Minimum Viable Model Blueprint

项目：Baryon-aided Time-dependent SIDM Accretion onto High-redshift Black-hole Seeds

日期：2026-07-04

## 1. 基线 SIDM 方程

沿用本地论文的 1D 球对称 SIDM conducting-fluid model：

```text
d rho_DM / dt + spherical advection = 0

d u_DM / dt + u_DM d u_DM / dr =
  - (1/rho_DM) d(rho_DM v_DM^2)/dr
  - G [M_DM(r,t) + M_BH(t)] / r^2

d(3 v_DM^2 / 2)/dt + advection =
  conduction
  - compressional heating/cooling
```

热传导使用 SMFP/LMFP 插值：

```text
kappa = 3/2 * [ sigma_m/(b v_DM)
               + 4 pi G/(a C rho_DM v_DM^3 sigma_m) ]^{-1}
```

黑洞 SIDM 吸积：

```text
dotM_DM = - 4 pi r_min^2 rho_DM(r_min) u_DM(r_min)
```

## 2. 第一版重子辅助：外加重子势场

最小修改是把重子质量分布加入动量方程：

```text
d u_DM / dt + u_DM d u_DM / dr =
  - (1/rho_DM) d(rho_DM v_DM^2)/dr
  - G [M_DM(r,t) + M_BH(t) + M_b(r,t)] / r^2
```

推荐的第一版重子 profile：

```text
Hernquist:
M_b(r,t) = M_b,tot(t) * r^2 / (r + a_b)^2
Phi_b(r,t) = - G M_b,tot(t) / (r + a_b)
```

时间增长模型：

```text
M_b,tot(t) = f_b M_halo [1 - exp(-t/t_b)]
```

初始测试参数：

```text
M_halo = 10^6, 10^7, 10^8, 10^9 M_sun
z = 10, 15, 20, 25, 30
M_seed = 1, 10, 100, 10^3 M_sun
f_b = 0, 0.01, 0.03, 0.05, 0.1, 0.16
a_b / r_s = 0.001, 0.003, 0.01, 0.03, 0.1
t_b = 1, 10, 30, 100 Myr
sigma/m = 1, 10, 30, 50, 100 cm^2/g
```

第一版的科学检验：

```text
Enhancement_DM(t) = dotM_DM(with baryons) / dotM_DM(no baryons)
Growth_DM(t) = M_BH(with baryons) / M_BH(no baryons)
t_reach(M_target) for M_target = 10^4, 10^5, 10^6 M_sun
```

## 3. 第二版重子辅助：重子吸积

黑洞质量增长拆成两项：

```text
dotM_BH = dotM_DM + dotM_b
```

重子吸积先用 Bondi-Eddington 截断：

```text
dotM_Bondi,b = 4 pi alpha G^2 M_BH^2 rho_b,inf
               / (c_s,b^2 + v_rel^2)^(3/2)

dotM_Edd = 4 pi G M_BH m_p / (epsilon_r sigma_T c)

dotM_b = min(dotM_Bondi,b, f_Edd dotM_Edd)
```

其中 `rho_b,inf` 和 `c_s,b` 第一版可由参数给定，第二版再从显式重子气体方程中求解。

关键输出：

```text
M_BH(t)
M_DM_accreted(t)
M_b_accreted(t)
f_DM,BH = M_DM_accreted / M_BH
Eddington-limited phase duration
dark-accretion-dominated phase onset time
```

## 4. 第三版：反馈参数化

反馈先用一参数模型，不立即做辐射流体：

```text
L_b = epsilon_r dotM_b c^2
dotE_fb = epsilon_f L_b
```

可选实现：

1. 加热中心重子气体，提高 `c_s,b`，降低 `dotM_Bondi,b`。
2. 让 baryon scale radius `a_b(t)` 随反馈能量增大。
3. 限制中心 baryon concentration，模拟 outflow。

推荐先做第 2 种：

```text
a_b(t) = a_b,0 [1 + E_fb(t) / E_bind,b]^eta
```

这能以最少自由度测试反馈是否抹掉重子增强效应。

## 5. 第一组数值实验

### Experiment 0：baseline reproduction

目的：复现本地论文中 NFW 与 SIS 的趋势。

输出：

- `M_BH(t)`
- `dotM_DM(t)`
- `rho_DM(r,t)`
- `u_DM(r,t)`
- `v_DM(r,t)`

### Experiment 1：static baryonic potential

目的：只打开静态 Hernquist baryon potential，看 SIDM 暗吸积增强多少。

对照：

```text
no baryon
f_b = 0.01, 0.05, 0.16
a_b/r_s = 0.001, 0.01, 0.1
```

### Experiment 2：growing baryonic potential

目的：检查重子形成时间是否重要。

对照：

```text
t_b = 1, 10, 30, 100 Myr
```

### Experiment 3：baryon accretion plus dark accretion

目的：连接 Feng/Yu/Zhong 的 baryon-aided dark Bondi scenario。

对照：

```text
f_Edd = 0.1, 1.0
epsilon_r = 0.1
M_seed = 1, 10, 100 M_sun
```

### Experiment 4：feedback robustness

目的：判断反馈能否终止重子辅助。

对照：

```text
epsilon_f = 0, 0.01, 0.05, 0.1
eta = 0.5, 1
```

## 6. 成功标准

第一篇工作的最低成功标准：

1. 能定量给出重子势场对 `dotM_DM` 的增强因子。
2. 能画出达到 `10^4, 10^5, 10^6 M_sun` 的参数空间。
3. 能说明哪些参数下增长主要来自暗物质，哪些参数下主要来自重子。
4. 能展示反馈是否显著收缩可行参数空间。

最强结果：

```text
在合理 high-z halo、baryon fraction、velocity-dependent SIDM cross section 下，
轻种子可以在 z ~ 4-11 前达到 LRD 所需质量，
且大部分最终质量来自暗物质吸积。
```

## 7. 下一步执行任务

1. 从本地论文中提取完整无量纲方程和参数标度，整理成代码需求。
2. 建立 Python 原型求解器目录结构。
3. 先实现静态重子势场模块，哪怕基线流体求解器尚未完整复现，也可以先验证 profile、质量、势场和单位换算。
4. 准备第一组 baseline 配置文件。

