# RIS 辅助 ISAC 仿真说明（公式详解版）

本文档对应当前 MATLAB 工程中的 RIS-ISAC 仿真流程，目标是给出与代码一致的数学建模、算法步骤和性能指标公式，便于论文写作与答辩讲解。

## 1. 系统与符号定义

### 1.1 维度与参数

- 基站天线数：M
- 通信用户数：K（每个用户单天线）
- RIS 单元数：N
- 总发射功率：P
- 雷达功率分配系数：eta，0 <= eta <= 1

功率划分：

$$
P_r = \eta P, \quad P_c = (1-\eta)P
$$

### 1.2 信道矩阵

- 直接通信信道：
  $$
  H_u \in \mathbb{C}^{K\times M}
  $$
- 直接雷达去程信道：
  $$
  h_{dt} \in \mathbb{C}^{M\times 1}
  $$
- BS 到 RIS 信道：
  $$
  G \in \mathbb{C}^{N\times M}
  $$
- RIS 到用户信道：
  $$
  H_{ru} \in \mathbb{C}^{K\times N}
  $$
- RIS 到目标信道：
  $$
  h_{rt} \in \mathbb{C}^{N\times 1}
  $$

### 1.3 设计变量

- RIS 相位向量：
  $$
  \phi=[\phi_1,\ldots,\phi_N]^T,\quad |\phi_n|=1
  $$
- RIS 对角矩阵：
  $$
  \Phi=\operatorname{diag}(\phi)
  $$
- 通信波束矩阵：
  $$
  W_c=[w_{c,1},\ldots,w_{c,K}] \in \mathbb{C}^{M\times K}
  $$
- 雷达波束向量：
  $$
  w_r \in \mathbb{C}^{M\times 1}
  $$

## 2. 复合信道与发射信号模型

### 2.1 用户等效信道

第 k 个用户的等效行向量来自“直达 + RIS 反射”：

$$
H_k = H_u + H_{ru}\Phi G \in \mathbb{C}^{K\times M}
$$

其中第 k 行记为 h_k（1 x M）。

### 2.2 雷达去程等效信道

$$
h_{\text{eff}}^T = h_{dt}^T + h_{rt}^T \Phi G \in \mathbb{C}^{1\times M}
$$

### 2.3 联合发射信号

设通信符号向量 s_c \in C^{Kx1}，雷达探测符号 s_r \in C，则发射信号为

$$
x = W_c s_c + w_r s_r
$$

总发射功率由 W_c 与 w_r 的缩放控制到 P。

## 3. 信道统计模型（与 generate_channels 一致）

### 3.1 路径损耗

统一形式：

$$
\beta(d)=10^{-3}d^{-\alpha}
$$

不同链路使用不同路径损耗指数 alpha_t, alpha_rt, alpha_k, alpha_rk, alpha_g。

### 3.2 Rician 模型

一般写法：

$$
H = \sqrt{\frac{\kappa}{1+\kappa}}H_{\text{LoS}} + \sqrt{\frac{1}{1+\kappa}}H_{\text{NLoS}}
$$

其中 H_NLoS 为复高斯随机项，kappa 为 Rician 因子。代码中 G、Hu、Hru 均采用该结构。

## 4. RIS 相位设计

## 4.1 固定目标相位（compute_phi: fixed_target）

代码使用：

$$
b = 2\,\operatorname{diag}(h_{rt}^T)\,\overline{G}\,h_{dt}
$$

并设置：

$$
\phi_n = e^{j\angle(b_n)},\quad n=1,\ldots,N
$$

即仅保留相位并满足单位模约束。其物理含义是对齐 RIS 反射分量与直达分量的相位。

### 4.2 启发式逐单元扫描（scan_phi）

对每个 RIS 单元 n 进行坐标下降式搜索：

1. 固定其余单元相位；
2. 在离散网格

$$
\mathcal{A}=\{0,\tfrac{2\pi}{12},\ldots,\tfrac{22\pi}{12}\}
$$

上尝试 cand_n=e^{ja}；
3. 每个候选相位下重新选择 eta、重算波束与目标函数；
4. 保留最优候选更新 phi_n。

目标函数可选：

- sumrate：最大化通信总速率
- snr：最大化雷达 SNR

## 5. 波束成形设计（design_w）

### 5.1 ZF 通信预编码

先构造伪逆形式：

$$
\widetilde{W}_c = H_k^H(H_kH_k^H)^{-1}
$$

代码写法等价于：

$$
\widetilde{W}_c = H_k'/(H_kH_k')
$$

然后逐列归一化（代码含数值稳定项 1e-12）：

$$
\widehat{w}_{c,k} = \frac{\widetilde{w}_{c,k}}{\|\widetilde{w}_{c,k}\|_2 + 10^{-12}}
$$

均分通信功率 Pc=(1-eta)P 给 K 个用户：

$$
w_{c,k} = \sqrt{\frac{(1-\eta)P}{K}}\,\widehat{w}_{c,k}
$$

### 5.2 MRT 雷达波束

基于去程等效信道 h_eff：

$$
w_r = \frac{\overline{h_{\text{eff}}}}{\|h_{\text{eff}}\|_2 + 10^{-12}}\sqrt{\eta P}
$$

这与最大比发射（MRT）一致。

## 6. 功率分配 eta 的可行性搜索（select_eta）

在离散集合

$$
\mathcal{E}=\{0,0.05,0.10,\ldots,1\}
$$

中从小到大搜索最小 eta，使雷达约束满足：

$$
g(\eta) \ge \gamma_t
$$

若 eta=1 仍不满足，则返回 eta=1 并给出不可行警告。

## 7. 通信性能公式（sum_rate）

第 k 个用户的 SINR：

$$
\operatorname{SINR}_k =
\frac{|h_k w_{c,k}|^2}
{\sum_{j\ne k}|h_k w_{c,j}|^2 + \sigma_k^2 + I_{r\to c,k}}
$$

其中雷达泄漏干扰项为可选：

$$
I_{r\to c,k} =
\begin{cases}
|h_k w_r|^2, & \text{若启用雷达对通信干扰建模} \\
0, & \text{否则}
\end{cases}
$$

总速率：

$$
R = \sum_{k=1}^{K}\log_2(1+\operatorname{SINR}_k)
$$

## 8. 雷达性能公式（radar_snr）

代码定义去程复合信道：

$$
h_{\text{forward}}^T = h_{dt}^T + h_{rt}^T\Phi G
$$

目标处雷达照射幅度：

$$
a_t = h_{\text{forward}}^T w_r
$$

单站假设下回程向量：

$$
h_{\text{backward}} = (h_{\text{forward}}^T)^T
$$

信号项：

$$
S = L\sigma_t^2|a_t|^2\|h_{\text{backward}}\|_2^2
$$

可选通信回波干扰项：

$$
I_{c\to r} =
\begin{cases}
L\sigma_t^2\sum_{k=1}^{K}|h_{\text{forward}}^T w_{c,k}|^2\|h_{\text{backward}}\|_2^2, & \text{启用} \\
0, & \text{否则}
\end{cases}
$$

最终雷达 SNR：

$$
g = \frac{S}{\sigma_r^2 + I_{c\to r}}
$$

## 9. 主流程（run_simple_isac）

对每个 N in {24, 32, ..., 64}：

1. 生成 Channel={hdt, Hu, hrt, G, Hru}
2. 计算四种 RIS 相位：fixed、random、scan-sumrate、scan-snr
3. 对每种相位：
   1. 计算 Hk=Hu+Hru*diag(phi)*G
   2. 用 select_eta 找到最小可行 eta
   3. 用 design_w 得到 (Wc, wr)
   4. 计算 R 与 g
4. 与无 RIS 基线比较

此外还绘制：

- 全局相位偏移 Delta 对性能的影响：phi = phi_base * e^{jDelta}
- 固定 N 下 eta-SNR 与 eta-Rate 权衡曲线

## 10. 关键结论（可直接用于答辩）

1. RIS 规模增大（N 上升）通常提升通信与感知两侧指标，体现无源阵列增益。
2. eta 增大时，雷达 SNR 单调增强趋势明显，而通信速率因功率让渡与同频泄漏而下降。
3. scan-sumrate 与 scan-snr 分别对应两种优化偏好，形成可解释的性能边界。
4. 该代码本质上实现了一个“离散 eta + 离散相位 + ZF/MRT 闭式波束”的低复杂度联合设计框架。

## 11. 可以直接写进论文的方法描述（模板）

可将本方法概括为如下优化问题：

$$
\max_{\Phi,\,W_c,\,w_r,\,\eta}\; R(\Phi,W_c,w_r,\eta)
$$

$$
	ext{s.t.}\quad g(\Phi,W_c,w_r,\eta)\ge\gamma_t,
\quad \sum_{k=1}^{K}\|w_{c,k}\|_2^2 + \|w_r\|_2^2 \le P,
\quad |\phi_n|=1,\,\forall n,
\quad 0\le\eta\le1
$$

工程实现上采用分层近似：

- 外层：逐单元离散扫描 RIS 相位；
- 中层：离散搜索最小可行 eta；
- 内层：给定 (phi, eta) 后使用 ZF + MRT 闭式解波束。

该分解避免了高维非凸联合优化的直接求解，复杂度可控且实现稳定。