---
layout: post
title: "Math Equations in Jekyll 2024"
date: 2025-03-01
tags: ['sasa','sasaas']
---

# Math Equations in Jekyll ✏️

考虑一个下行多用户MISO系统，基站配有 $ M $ 个天线，服务 $ N $ 个单天线用户。用户 $K$ 的符号为 $x_k$，给其的precoding为 $ \mathbf{w}_k \in \mathbb{C}^{M \times 1} $。基站发出的信号为

$$
\mathbf{s}=\sum_{i=1}^N \mathbf{w}_i x_i  \tag{1}
$$

用户$k$ 的接收信号为

$$
\mathbf{y}_k = \mathbf{h}_k^H \mathbf{s}  \tag{2}
$$

基站到用户$k$的信道为 $\mathbf{w}_k \in \mathcal{C}^{M\times 1}$ ， 用户 $k$ 的SINR 为

$$
\gamma_{k}=\frac{\left|\mathbf{h}_{k}^{H} \mathbf{w}_{k} \right|^2}{\sum_{j \neq k} \left|\mathbf{h}_{k}^{H}\mathbf{w}_{j} \right|^2 + \sigma^2},
$$

**优化问题**

目标函数为最小化基站总功率，约束条件为每个用户的SINR约束，

$$
\begin{equation}

\begin{aligned}

    & \text{P1}: \,\, \min_{\{\mathbf{w}_{k}\}\forall k} \sum_{k=1}^U \|\mathbf{w}_{k}\|_2^2, \quad \text{s.t.} \\
    & \text{C1}:\quad \frac{\left|\mathbf{h}_{k}^{H} \mathbf{w}_{k} \right|^2}{\sum_{ j \neq k} \left|\mathbf{h}_{k}^{H}\mathbf{w}_{j} \right|^2 + \sigma^2} \geq \gamma_k.
    
\end{aligned}

\end{equation}
$$


此问题可以直接由cvx求解

你可以写 **行内数学公式**：

   000  Euler’s formula: $ [1] \equiv [3] $ : It follows from the monotone gradient condition for convexity of $ g(x) $ , i.e., $ g(x) $  is convex if and only if $ (\nabla g(x) - \nabla g(y))^T(x-y) \ge 0,~\forall x,y. $  SINR is $ \gamma_{s,k}=\mathbf{h}{s,k}^{H} $ 

```
N = 8;                                                                      % Number of users
M = 5;                                                                      % Number of antennas at each satellite
K=M;
lambda = 0.1;                                                                 % Weighting factor for alpha penalty
mu=0.01;
P = 10;                                                                     % Power constraint
N0_dBm_Hz = -173;                                                           % Noise power spectral density in dBm/Hz
B = 1e6;                                                                    % Bandwidth in Hz (adjust if necessary)
N0_W_Hz = 10^((N0_dBm_Hz - 30) / 10);                                       % Convert N0 from dBm/Hz to Watts/Hz
sigma2 = 1;                                                       % noise power
```


你也可以写 **块级数学公式**：
$$
E = mc^2
$$

$$
\gamma_{s,k}=\frac{\left|\mathbf{h}_{s,k}^{H} \mathbf{v}_{s,k} \right|^2}{\sum_{j\in \mathcal{S}_s, j \neq k} \left|\mathbf{h}_{s,k}^{H}\mathbf{v}_{s,j} \right|^2 + \sigma^2},
$$
### 📌 **更多公式示例**
1. 级数求和：
   $$
   S_n = \sum_{i=1}^{n} a_i
   $$
2. 定积分：
   $$
   \int_{a}^{b} x^2 \,dx = \frac{b^3}{3} - \frac{a^3}{3}
   $$

3. 矩阵：
   $$
   A = \begin{bmatrix} a & b \\ c & d \end{bmatrix}
   $$

4. 逻辑符号：
   $$
   P(A \cap B) = P(A) P(B \mid A)
   $$
