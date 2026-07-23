import math

class Binom:
    """
    二项分布

    Property
    ------
    n : int
        重复试验的次数
    p : float
        每次试验成功的概率，必须满足 0 <= p <= 1

    Method
    ------
    pmf(k)
        计算概率质量函数 P(X = k)
    cdf(k)
        计算累积分布函数 P(X <= k)
    expectation()
        返回期望 E[X]】
    variance()
        返回方差 Var[X]
    std_dev()
        返回标准差 Std[X]
    """

    def __init__(self, n, p):
        """
        初始化二项分布对象

        Parameters
        ------
        n : int
            总的独立重复试验次数
        p : float
            每次成功的概率,0 <= p <= 1
        """
        if not isinstance(n, int) or n <= 0:
            raise ValueError('试验次数n必须为正整数')
        if not (0 <= p <=1):
            raise ValueError('成功概率必须在范围[0,1]内')

        self.n = n
        self.p = p

        self.E = self.expectation()
        self.Var = self.variance()
        self.Std = self.std_dev()


    def pmf(self, k:int) -> float:
        """
        概率质量函数（PMF）

        参数
        ------
        k : int
            成功的次数

        返回
        ------
        float : P(X = k) 的概率
        """
        if k < 0 or k > self.n:
            return 0.0

        comb = math.comb(self.n, k)
        prob = comb * (self.p ** k) * ((1 - self.p) ** (self.n - k))
        return prob


    def cdf(self, k:int) -> float:
        """
        累积分布函数（CDF）

        参数
        ------
        k : int
            成功次数阈值

        返回
        ------
        float : P(X <= k) 的累积概率
        """
        if k < 0:
            return 0.0
        if k > self.n:
            return 1.0

        total = 0.0
        for i in range(0, k+1):
            total += self.pmf(i)
        return total


    def expectation(self) -> float:
        """
        返回期望值 E[X] = n * p
        """
        return self.n * self.p

    def variance(self) -> float:
        """
        返回方差 Var[X] = n * p * (1 - p)
        """
        return self.n * self.p * (1 - self.p)

    def std_dev(self) -> float:
        """
        返回标准差 Std[X] = sqrt(Var[X])
        """
        return math.sqrt(self.variance())