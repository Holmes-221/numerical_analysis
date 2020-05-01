***验证Picard迭代，Newton迭代***

***要求：***

***构造函数以及迭代方法，判定迭代收敛速度优略***

***找到Newton迭代法与Picard迭代法的优势***

***找到Newton迭代法的缺陷***

**构造函数**$\mathbf{f}\left( \mathbf{x} \right)\mathbf{=}\mathbf{x}^{\mathbf{3}}\mathbf{- x - 1}$**，迭代初始值为x0=1.5**

**<span style="font-variant:small-caps;">Newton迭代法</span>**

**1.收敛状况**

Newton迭代格式为 []{#_Hlk37433039
.anchor}$x_{n + 1} = \varphi\left( x_{n} \right) = x_{n} - \frac{f\left( x_{n} \right)}{f'\left( x_{n} \right)}$，当$\left| \varphi\left( x \right) \right| \leq L < 1$时，Newton迭代法收敛，L越趋于0，收敛效果越好；L越趋于1，收敛效果越差；当L&gt;1时，Newton迭代法发散。按照原方程的解根的重数，一般分为两种情形，单根和重根。

***（1）单根***

迭代格式为$x_{n + 1} = \varphi\left( x_{n} \right) = x_{n} - \frac{x_{n}^{3} - x_{n} - 1}{3x_{n}^{2} - 1}$，由于有

$|\varphi'\left( x0 \right)| = |\frac{f\left( x0 \right)f''\left( x0 \right)}{f'\left( x0 \right)^{2}}| = 0.2381852551984877 < 1$，故Newton迭代法收敛，为了方便计算，取L=0.3。当迭代6次时，精度达到$10^{- 4}$，其值为1.32471796，在单根的情况下，由于$\varphi\left( x* \right) \neq 0,\varphi^{'}\left( x* \right) = 0,\varphi''(x*) \neq 0$，故Newton迭代格式二阶收敛

***（2）重根***

构造函数$f = \left( x - 1 \right)\left( x - 2 \right)^{2}$，其Newton迭代格式为
$x_{n + 1} = \varphi\left( x_{n} \right) = x_{n} - \frac{x_{n}^{2} - 3x_{n} + 2}{3x_{n}}$，取初值x0=2.5，$\left| \varphi^{'}(x0) \right| = \left| \frac{f\left( x0 \right)f^{''}(x0)}{f^{'\left( x0 \right)^{2}}} \right| = 0.6122448979591837 < 1$，取L=0.7。故Newton迭代格式收敛，由于$\varphi\left( x* \right) \neq 0,\varphi^{'}\left( x* \right) = \frac{1}{2} \neq 0$，因此Newton迭代法一阶收敛。

**2.优势**

Newton迭代法是Picard迭代法的一种特殊形式，不需要另外构造迭代格式，由于Picard迭代函数的未知性，很难找到收敛的迭代序列，更别说二阶收敛乃至更高阶的收敛序列。

在单根的情况下，迭代格式至少是二阶收敛，收敛速度较快。

**3.缺陷**

Newton迭代法是局部收敛，对初值的要求高，如果初值选择不当可能导致不收敛的情况。例如在上述单根情况中$\varphi\left( x_{n} \right) = x_{n} - \frac{x_{n}^{3} - x_{n} - 1}{3x_{n}^{2} - 1}$，当初值x0=-1时，$|\varphi'(x0)| > 1$，故迭代格式不收敛。而当初值选为1.5时，迭代格式收敛。

需要事先判断方程的解是否为重根，因为在重根情况下，迭代格式只是一阶线性收敛，收敛速度较慢。

**4.改进**

如果是重根，可以改进Newton迭代格式$x_{n + 1} = \varphi\left( x_{n} \right) = x_{n} - m\frac{f\left( x_{n} \right)}{f'\left( x_{n} \right)}$，其中m为根的重数。

此时$\varphi\left( x* \right) \neq 0,\varphi^{'}\left( x* \right) = 0$，此时迭代格式至少二阶收敛。

**<span style="font-variant:small-caps;">Picard迭代法</span>**

**1.收敛状况**

Picard迭代格式为
$x_{n + 1} = \varphi\left( x_{n} \right)$，并且对任意x1，x2在初始值的领域内，存在$0 < L < 1,$都有$\left| \varphi\left( x1 \right) - \varphi\left( x2 \right) \right| < L\left| x1 - x2 \right|$成立。其中L越趋于0，收敛效果越好；L越趋于1，收敛效果越差；当L&gt;1时，迭代格式发散。从这里可知Newton迭代是一种特殊的Picard迭代格式。特别地，L可以用$\left| \varphi^{'}\left( x \right) \right|$替代。

构造Picard迭代序列 $\varphi 1\left( x_{n} \right) = x_{n}^{3} - 1$ ,
$\varphi 2\left( x_{n} \right) = \left( x_{n} + 1 \right)^{\frac{1}{3}}$；初始条件x0=1.5

先讨论第一种迭代序列$\varphi 1\left( x_{n} \right) = x_{n}^{3} - 1$，$\varphi 1^{'}\left( x \right) = 3x^{2}$，易知${|\varphi}^{'}\left( x0 \right)| = \frac{27}{4} > 1$，故迭代序列发散。对于第二种迭代序列$\varphi 2\left( x_{n} \right) = \left( x + 1 \right)^{\frac{1}{3}}$，有$\varphi 2^{'}\left( x \right) = \frac{1}{3}\left( x + 1 \right)^{- \frac{2}{3}}$，当x0=1.5时，$\left| \varphi 2^{'}\left( x0 \right) \right| = 0.18096117443966045 < 1$，因此迭代序列收敛。为方便计算，不妨记L=0.2；达到精度$\ 10^{- 4}$
仅需要4次Picard迭代，值为1.324760011。

**2.优势**

不需要讨论根的分布情况，由于迭代格式是未知的，有可能出现收敛三阶收敛甚至更高阶收敛的情况。从上面的例子我们可以看出，对于求解同一个方程$f\left( x \right) = x^{3} - x - 1$，x0=1.5，Newton迭代格式$\varphi\left( x_{n} \right) = x_{n} - \frac{x_{n}^{3} - x_{n} - 1}{3x_{n}^{2} - 1}$和Picard迭代序列$\varphi 2\left( x_{n} \right) = \left( x_{n} + 1 \right)^{\frac{1}{3}}$，在达到相同精度($10^{- 4}$)的前提下，Picard迭代法只需要迭代4次，而Newton迭代法需要迭代6次。
