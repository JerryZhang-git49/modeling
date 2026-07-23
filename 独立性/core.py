import math
from functools import reduce  # reduce函数用于将一个函数连续作用在序列的元素上


def series(num,p):
    """
        用于一键计算串联系统的算放框架————必须全部连通
    Parameters
    ----------
    num: int
        该系统中的元素数量
    p: float | list[float]
        故障率，若为一个[0, 1]之间的常数，表示各个零部件故障率相同

        若为一个列表，则列表长度应该与num的值相同，表示每个零件对应的故障率

    Returns
    -------
    out: float
        该系统的可靠性

    Examples
    --------
    >>> series(3, 0.1) # 三个零件，每个零件的故障率都是0.1
    0.729

    >>> series(3, [0.1, 0.2, 0.3])  # 三个零件，故障率分别为0.1、0.2、0，3
    0.504
    """
    if isinstance(p,float):
        if p < 0 or p > 1:
            raise ValueError('p must be beyween 0 and 1')
        p = [p] * num   # 如果是单个浮点数，就复制出num个相同的值
    if len(p) != num:
        raise IndexError('length of p must be equal to num')

    p = [1 - i for i in p]   # 将故障率转化为可靠率
    return round(reduce(lambda x,y:x*y, p), 5)  # 将可靠率连续相乘得到串联下的系统可靠率

def parallel(num, p):
    """
        用于一键计算并联系统的可靠性————有一个连通即可
    Parameters
    ----------
    num: int
        该系统中的元素数量
    p: float | list[float]
        故障率，若为一个[0, 1]之间的常数，表示各个零部件故障率相同

        若为一个列表，则列表长度应该与num的值相同，表示每个零件对应的故障率

    Returns
    -------
    out: float
        该系统的可靠性

    Examples
    --------
    >>> parallel(3, 0.9) # 必须全部坏，该系统才会坏，所以这里把故障率设为0.9，那么全坏的概率就是0.9^3=0.729，对应可靠性0.271
    0.271

    >>> parallel(3, [0.9, 0.8, 0.7])  # 可以看到，概率与串联系统取反的情况下，所得到的值也是1-串联系统的概率
    0.496
    """
    if isinstance(p, float):
        if p < 0 or p > 1:
            raise ValueError('p must be between 0 and 1')
        p = [p] * num
    if len(p) != num:
        raise IndexError('length of p must be equal to num')

    return round(1 - reduce(lambda x,y:x*y, p), 5)   # 1 - 故障率连续相乘


def vote(num, p, r):
    """
            用于一键计算表决系统的可靠性————有r个以上或比例大于r连通时
        Parameters
        ----------
        num: int
            该系统中的元素数量

        p: float
            故障率，为一个[0, 1]之间的常数，表示各个零部件故障率相同

        r: int| float
            若为整数类型，则表示该系统中最少需要多少台机器正常才正常，此时r<=num

            若为小数类型，则表示该系统中最少需要多少占比的机器正常才正常，算得结果向上取整，此时 r \in [0, 1]

        Returns
        -------
        out: float
            该系统的可靠性

        Examples
        --------
        >>> vote(3, 0.1, 1)
        0.999

        >>> vote(3, 0.1, 2)
        0.972

        """
    if isinstance(p, float):
        if p < 0 or p > 1:
            raise ValueError('p must be between 0 and 1')
        p = [p] * num
    else:
        raise ValueError('p必须是一个数字')

    def binom(k, n, p2):
        # 根据k,n,p计算二项分布概率值
        return math.comb(n, k) * (p2**k) * (1-p2)**(n-k)

    result = 0

    # 若r为占比，则向上取整求数量
    if isinstance(r, float) and 0<r<1:
        r = math.ceil(num*r)

    # 循环求概率（筛选连通数量大于r的）
    # 尽可能减少循环次数
    if r > num/2:
        # 如果需要的正常设备占多数，则直接计算从r到n所有正常情况的概率和
        for i, _p in zip(range(r, num + 1), p[r-1:]):
            result += binom(i, num, 1-_p)
    else:
        # 否则先计算失败的概率（低于r个正常），再用 1 减去得到可靠性
        for i, _p in zip(range(r), p[:r]):
            result += binom(i, num, 1-_p)
        result = 1 - result

    return round(result, 5)