#### 实现方程组求解中的高斯消去法、列主元高斯消去法、LU分解法、列主元LU分解法要求：

#### 构造方程组，实现高斯消去法、列主元高斯消去法、LU分解法、列主元LU分解法

#### 发现实现高斯消去法、列主元高斯消去法区别（为什么要用列主元高斯消去法）

#### 发现实现LU分解法、列主元LU分解法区别（为什么要用列主元LU分解法）构造方程组 $\begin{pmatrix}
\mathbf{1} & \mathbf{1} & \mathbf{1} \\
\mathbf{1} & \mathbf{3} & \mathbf{- 2} \\
\mathbf{2} & \mathbf{- 2} & \mathbf{1} \\
\end{pmatrix}\begin{pmatrix}
\mathbf{x1} \\
\mathbf{x2} \\
\mathbf{x3} \\
\end{pmatrix}\mathbf{=}\begin{pmatrix}
\mathbf{6} \\
\mathbf{1} \\
\mathbf{1} \\
\end{pmatrix}$\
为了方便计算，以下记系数矩阵$\mathbf{A =}\begin{pmatrix}
\mathbf{1} & \mathbf{1} & \mathbf{1} \\
\mathbf{1} & \mathbf{3} & \mathbf{- 2} \\
\mathbf{2} & \mathbf{- 2} & \mathbf{1} \\
\end{pmatrix}\mathbf{,}$

#### $\mathbf{}\mathbf{X}\mathbf{=}\begin{pmatrix}
\mathbf{x1} \\
\mathbf{x2} \\
\mathbf{x3} \\
\end{pmatrix}\mathbf{,}\mathbf{b =}\begin{pmatrix}
\mathbf{6} \\
\mathbf{1} \\
\mathbf{1} \\
\end{pmatrix}$ 

线性方程组的直接解法

线性方程组可表示为Ax=b，如果det(A) $\neq$
0，在此情形下x有唯一解。在线性代数中我们学过一种解线性方程的 Cramer's
Rule，在理论上它是一种可行的方法，但在实际计算中工作量巨大，因此一般不不采用这种方法求解方程组。目前在计算机上常使用的线性方程组的数值解法一般有两种：直接法和迭代法。直接法一般用于求解稠密矩阵中的低阶方程组。迭代法一般用于求解稀疏矩阵的高阶方程组。

##### Gauss消去法

如果系数矩阵A可逆，可将A经过初等行变换线性变换成上三角矩阵U(行阶梯形矩阵)，该方程组与原方程组同解。Gauss消去法可分为两个步骤：1.消元过程
2.回代过程。其中U的对角线元素$u_{\text{ii}}$（主元素）不等于0（A的对角线元素全不为0）。消元过程的算法计算量为$\frac{n\left( n - 1 \right)\left( 2n + 5 \right)}{6}$，回代过程中算法的计算量约为$\frac{n\left( n + 1 \right)}{2}$。因此Gauss消去法的计算量为$O\left( n^{3} \right),$且主要计算量体现在消元过程。从最后一个方程向上回代则可以求解出
$x_{1},x_{2},x_{3},\ldots,x_{n}$。另外Gauss消去法的充要条是系数矩阵A的所有顺序主子式不为0。
例如在本题中，利用Gauss消去法可得到：$(A|b) \rightarrow \begin{pmatrix}
1 & 1 & 1 & 6 \\
0 & 2 & - 3 & - 5 \\
0 & 0 & - 7 & - 21 \\
\end{pmatrix}$
利用第一行消去第一列的元素；利用第二行消去第二列的元素………以此类推。再从最后一个方程回代，可解出
$x_{1} = 1,x_{2} = 2,x_{3} = 3$。

特别地，如果方程组的系数矩阵A是严格对角占优(对角线元素的绝对值大于该行其他元素值之和)或是对称正定的，则可采用Gauss消去法，并且算法稳定。

##### 列主元Gauss消去法

Gauss消去法虽然能求解出一些线性方程，担任由于一些不可避免的缺陷。

1.如果在化成行阶梯形矩阵时出现主对角线元素为0(存在顺序主子式为0)，那么Gauss消去法失效。

2.如果主对角线元素$u_{\text{ii}}$特别小时（趋于0），由于它是作为除数，由误差分析可以得知此时会造成非常大的误差(舍入误差的存在不能保证计算的过程仍是数值稳定的)。例如
$\begin{pmatrix}
0.03 & 58.9 \\
5.31 & - 6.10 \\
\end{pmatrix}\begin{pmatrix}
x_{1} \\
x_{2} \\
\end{pmatrix} = \begin{pmatrix}
59.2 \\
47.0 \\
\end{pmatrix}$
使用3位浮点数计算。易知方程组的准确解为$x_{1} = 10,x_{2} = 1$。
如果采用Gauss消去法求解，会得到数值解$x_{1} = - 10.0,x_{2} = 1.01$，显然这个解的严重偏离了精确解，并不是我们想要的结果。

由此可见对一些系数矩阵A无法采用Gauss消去法或者用Gauss消去法会导致误差严重增大甚至出现错误。

如果我们采用列主元Gauss消去法，则避免了Giass消去法的局限性。它是在Gauss消去的过程中做了一点改进，即在对矩阵的主元$a_{\text{kk}}^{\left( k \right)}$进行除法前选择剩余行中第k列元素绝对值最大的值$a_{\text{jk}}^{\left( k - 1 \right)}$,然后第k行与第j行互换，一直重复此操作。与Gauss消去法相比较，列主元Gauss消去法是对系数矩阵A做一系列的初等行变换
(左乘一个可逆矩阵P)。再次题中$\begin{pmatrix}
A| & b \\
\end{pmatrix} = \begin{pmatrix}
2 & - 2 & 1 & 1 \\
0 & 4 & - \frac{5}{2} & \frac{1}{2} \\
0 & 0 & \frac{7}{4} & \frac{21}{4} \\
\end{pmatrix}$，解得$X = \begin{pmatrix}
3 \\
2 \\
1 \\
\end{pmatrix}$。从上边可以看出，列主元消去法避免了顺序消去法解不稳定的问题。

##### LU分解法 

如果A的各阶顺序主子式不为0，在Gauss消元法中我们知道左乘一个可逆矩阵可以将系数矩阵A化为上三角矩阵U，很容易证明那可逆矩阵是一个单位下三角矩阵L，即有A
= LU，因此原方程可化为
LUX=b(Doolittle分解法)。记$UX = y,Ly = b$，则可解出未知量X。同样，还可以将A分解成为下三角和单位上三角矩阵(Crout分解法)，本题是采用的Doolittle分解法。

$A = \begin{pmatrix}
1 & 1 & 1 \\
1 & 3 & - 2 \\
2 & - 2 & 1 \\
\end{pmatrix} = \begin{pmatrix}
1 & 0 & 0 \\
1 & 1 & 0 \\
2 & - 2 & 1 \\
\end{pmatrix}\begin{pmatrix}
1 & 1 & 1 \\
0 & 2 & - 3 \\
0 & 0 & - 7 \\
\end{pmatrix} = LU$。

其紧凑格式为$\begin{pmatrix}
1 & 1 & 1 \\
1 & 2 & - 3 \\
2 & - 2 & - 7 \\
\end{pmatrix}$ $\begin{pmatrix}
1 & 0 & 0 \\
1 & 1 & 0 \\
2 & - 2 & 1 \\
\end{pmatrix}y = \begin{pmatrix}
1 \\
1 \\
6 \\
\end{pmatrix},\begin{pmatrix}
1 & 1 & 1 \\
0 & 2 & - 3 \\
0 & 0 & - 7 \\
\end{pmatrix}X = y$，解出$x_{1} = 3,x_{2} = 2,x_{3} = 1$。
与Gauss消去法相比，他能解所有相同系数矩阵的方程组，而Gauss消去法只能解只一个方程组，当遇见相同系数矩阵的方程组求解时效率比较低下。

##### 列主元的LU分解法

从列主元Gauss消去法我们可以得知，相当于对系数矩阵做一个初等行变换(左乘一个可逆矩阵P)。这么做的原因和列主元Gauss消去法的原因大致相同。1.在化成行阶梯形矩阵时出现主对角线元素为0，那么Gauss消去法失效；2.如果主对角线元素\$a\_{ii}特别小时（趋于0），而被除数又特别大时。由于它是作为除数，由误差分析可以得知此时会造成非常大的误差(在做除法时应避免绝对值很小的数作为除数)。求解方法和列主元Gauss消去法也大致相同。

$AX = b \rightarrow PAX = Pb \rightarrow LUX = Pb$，记$UX = y,Ly = Pb$。
然后从最后一个方程回代即可解出X。
