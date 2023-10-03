## Notes For "Inverse Problems for Partial Differential Equations"

### Ch1 Inverse Problems

​	本章主要介绍一些反问题，包括

1. 反引力问题（及相似的反位势问题）——由势函数的梯度函数确定其“质量”分布
2. 反电传导（Conductivity）问题——与电位势相关
3. 反散射问题——通过远场（Far Field）确定散射物体 Scatter
4. 扫描断层技术的数学表示——通过在各个方向超平面上的积分确定整体的分布，如著名的Radon变换
5. 反谱（Spectral）问题——给定方程，定义域会决定微分算子的谱，我们要做的就是从谱反推出定义域的形状

### Ch2 不适定问题与正则化 Regularization

​	给一个正问题算子
$$
A: \mathcal{X} \to \mathcal{Y}, Ax=y
$$
我们想要找一个逆的算子
$$
A^{-1}: A\mathcal{X}\to \mathcal{X}, Ax\to x
$$
一般而言，这是做不到的.正则化的思路就是找一族适定（连续）的算子
$$
R_\alpha :A\mathcal{X}\to \mathcal{X}, \lim_{\alpha \to 0}R_\alpha Ax=x
$$
我们称这族算子为原来算子的正则子

##### 一大类正则子的构造

​	可以通过增加一个正则项 $\alpha \mathcal{M}$ 来构造.其中 $\mathcal{M}$ 是一个非负下半连续的函数，可以选择一范数、二范数、香农熵等等

##### 正则化参数 $\alpha$ 的选取

另一个比较有意思的点是我们虽然可以通过使 $\alpha\to 0$ 取极限得到 $A^{-1}$ 在一点的取值，但是实际上由于测量误差，舍入误差的存在，这个点的位置也不是精确的.所以事实上我们选取一个不大不小的 $\alpha$ 是最好的.

### Ch3 柯西问题的唯一性与稳定性

主要分析不同方程的柯西问题的反演问题

#### 3.1 抛物方程

​	一般形式：$\partial_t u + A u =0\ on\ \Omega \times (0.T)$

###### 特征函数法

​	对 $A=-\it{div}(a\nabla)+c$ 来说，其逆算子是自共轭的紧算子，在 $H^2\cap H_0^1(\Omega)$ 中有一组由其特征向量组成的基. 由此可以表示出唯一性，进一步，借助名为 Quasi-reversibility 的操作，我们给圆方程加上一个关于 A 的高阶项，得到一个对特征函数方法稳定的 regularizer  .

###### 对数凸性法 The method of the logarithmic convexity

​	记 $f(t)=\|u(t)\|_2^2(\Omega)$ ，而 $F(t)=ln f(t)$ . 对很多算子 $A$，例如 $A=-\it{div}(a\nabla)+c$ ，满足 $F''\geq 0$ ，即 $F$ 具有凸性. 由此可以利用 $F(0)$ 与 $F(T)$ 限制 $F(t), t\in (0,T)$ 的值.

​	该方法可以扩展到不严格的方程 $\|\partial_t u + A u\|\leq \alpha \|u\|$.

​	当 $A$ 是希尔伯特空间(基域为复数域)中的线性算子，且可以分解为对称算子和反对称算子的和时， $A=A_++A_-$ ，若满足一定的正则条件，对 $\|A_-u\|^2$ 以及 $\partial_t (A_+ u, u)=\partial_t (A u, u)$ 有一个合适的限制，我们可以得到一个很漂亮的结果
$$
\|u(t)\|\leq C_1 \|u(0)\|^{(1-\lambda)} \|u(T)\|^{\lambda}
$$
其中 $C_1,\lambda$ 是由 $\alpha$ 决定的常数，当 $\alpha=0$ 时，可以取 $C_1=1, \lambda =t/T$.

​	**Example** The method of logarithmic convexity 可以用在有界区域上定义的 
$$
Au= -\sum\partial_k(a_{jk} \partial_j u) + \sum b_j \partial_j u + cu
$$
其定义域为 $D(t)= H_0^1(\Omega)\cap H^2(\Omega)$. $a_{jk}, \partial_t a_{jk}, b, {\rm div} b \in C^1(\overline{\Omega} \times [0,T])$， $c, \partial_t c \in L_\infty (\Omega \times [0, T])$.

​	我们将 $A$ 分解为
$$
A_+ u = -\sum\partial_k(a_{jk} \partial_j u) + (-\frac{1}{2} {\rm div\ } b + c) u
$$

$$
A_- u = \sum \partial_j(b_j u) - \frac{1}{2} ({\rm div}\ b)u
$$

验证其符合定理中的条件，应用定理，对方程 $\partial_t u + Au=0$，其解应该满足条件
$$
\|u(t)\|\leq \|u(0)\|^{(1-\frac{t}{T})} \|u(T)\|^{\frac{t}{T}}
$$

###### Semigroups

完全不懂。。

#### 3.2 Carleman Estimates

​	在本书中， Carleman 估计提供了一个 $\|Au\ w\|_2$ 的下界，即 $Au$ 的带权重2范数有下界，可以借此得到一定条件下的唯一性. 在 Lerner 的 Carleman Inequalities: An Intro and More 中， Carleman 估计则是作为解延拓的一个手段.(也可以用来做唯一性)

​	Lerner 提到，在 Carleman 估计之前，方程解的唯一性结果要么依赖于方程的双曲性，要么借助 Cauchy-Kovalevskaya 定理、Holmgren 定理（对微分算子有很强的解析性要求 strong analyticity structure）.

###### 一般的微分算子的情形

​	给定微分算子$A(x,\partial)=\sum_{|\alpha:m|\leq 1} a_\alpha \partial^\alpha$，我们要找一个满足所谓 strong pseudoconvexity 条件的 $C^2$ 函数 $\phi$,定义一个权重函数 $w(x)=\exp(\tau \phi(x))$. 那么存在一个常数 $C$ 使
$$
\tau \int_\Omega |\partial^\alpha u|^2 w^2  dx\leq C\int_\Omega |Au|^2 w^2 dx, |\alpha:m|<1,
$$
即 $\tau \|\partial^\alpha u\ w(x;\tau)\|_2^2\leq C \| Au\ w(x;\tau)\|_2^2$对所有 $\tau>C$ 成立. (对右侧使用 Holder 不等式放缩，除至左边，可以得到 $\|Au\|_2$ 的下界估计，可以用来导出矛盾，得到唯一性).

###### 二阶微分算子的情形

​	二阶时，可以将对 $\phi$ 的要求降低.

Theorem : 设 $A$ 是一个二阶微分算子，其二阶部分 $A_m(x;\xi)$ 的系数为实值函数. 区域是非特征的，即 $A(x;\nabla \psi(x))\neq 0, x\in \overline{\Omega}$.

​	对函数 $\psi\in C^2(\overline{\Omega}), \partial \Omega \in C^2$，若满足
$$
A_m(x;\xi)=0,\sum(\partial A_m/\partial \xi_j)\partial_j \psi =0, \xi\neq 0,
$$
能推出
$$
\sum (\partial_j \partial_k\psi)\frac{\partial A_m}{\partial\xi_j} \frac{\partial A_m}{\partial\xi_k}\\
+\sum \partial_j\psi (\partial_k \frac{\partial A_m}{\partial\xi_j}\frac{\partial A_m}{\partial\xi_k} -\partial_k A_m \frac{\partial^2 A_m}{\partial\xi_j\partial\xi_k})>0,
$$
其中 $\partial_j,\partial_k$ 代表对 $x$ 的坐标求导. 引入函数
$$
\phi= e^{\sigma \psi}.
$$
那么存在常数 $C_1(\sigma),C_2$ 使
$$
\tau^{3-2|\alpha|} \int_\Omega |\partial^\alpha u|^2 w^2  \leq C_1\left(
\int_\Omega |Au|^2 w^2 +\int_{\partial \Omega} (\tau|\nabla u|^2+\tau^3 |u|^2)w^2
\right).
$$
对所有 $C_2<\sigma, C_1<\tau,|\alpha|\leq 1$ 以及 $u\in H^2(\Omega)$ 成立.

##### 由 Carelman 估计得到唯一性与稳定性

直接应用 Carelman 估计，可以得到估计
$$
\|\partial^\alpha u\|_2(\Omega_\epsilon) \leq C(F+M^{1-\kappa}F^{\kappa})\ when\ |\alpha:m|<1,
$$
其中 $F=\|f\|_2(\Omega) +\sum \|g_j\|_{(m_1-j-1/2)}(\Gamma)$， $M$ 是 $\|\partial^\alpha u\|_2(\Omega)$ 关于 $|\alpha:m|<1$ 求和.

#### 3.3 由 Carleman 估计具体计算一些椭圆方程与抛物方程

















