"""
    波利亚模型通解

b : int
    抽样类型1的数量
r : int
    抽样类型2的数量
m : int
    抽到m个类型1
n : int
    抽到n个类型2
c : int
    每次抽完之后,再放入c个同类型的样本
d : int
    每次抽完之后,再放入d个不同类型的样本
"""

from math import comb

def sampling_without_replacement(b, r, m, n, /) -> float | None:
    """
    无放回抽样(超几何分布)
    只涉及两种类型的抽样计算，不涉及多维超几何分布的计算  \n
    例如,检测次品率,在100个零件中进行抽取,其中好零件有90个,坏零件有10个,
    现在抽取了10次,抽取的结果分别是8个好零件和2个坏零件 \n
    那么计算的结果就是sampling_without_replacement(90, 10, 8, 2)
    Parameters
    ----------
    b : int
        抽样类型1的数量
    r : int
        抽样类型2的数量
    m : int
        抽到m个类型1, m < b
    n : int
        抽到n个类型2, n < r
    """
    # 参数校验
    if not all(isinstance(x,int) and x >= 0 for x in [b,r,m,n]):
        raise ValueError('所有参数必须是非负整数')

    if m > b or n > r:
        raise ValueError(f'抽样数不能超过总体:m={m} > b={b} 或 n={n} > r={r}')

    if m+n > b+r:
        raise ValueError(f"抽取总数 {m+n} 不能超过总体 {b+r}")

    # 计算概率
    return comb(b,m) * comb(r,n) / comb(b+r,m+n)


def sampling_with_replacement(b:int, r:int, m:int, n:int, /) -> float:
    """
    放回抽样(二项分布)
    从包含b个类型1和r个类型2的总体中,有放回地抽取(m+n)次，
    恰好得到m个类型1和n个类型2的概率。
    Parameters
    ----------
    b : int
        抽样类型1的数量
    r : int
        抽样类型2的数量
    m : int
        抽到m个类型1
    n : int
        抽到n个类型2
    """
    # 参数校验
    if not all(isinstance(x,int) and x >= 0 for x in [b.r.m.n]):
        raise ValueError('所有参数必须是非负数')

    if m+n == 0:
        return 1.0   # 抽0次，必然得到0个

    p = b / (b+r)
    return comb(m+n,m) * (p ** m) * ((1-p) ** n)


def infection(b:int, r:int, m:int, n:int, c:int, /) -> float:
    """
    infection(b, r, m, n, c, /) -> float \n
    计算方法：C_(m+n)^m\prod{\prod^{m-1}{i=0}(r+ic)\prod^{n-1}{j=0}(b+ic)}{\prod^{m+n-1}{k=0}(b+r+kc)} \n
    波利亚模型中的传染病模型，每次发现一个传染病都会增加再感染的概率，每次发现一个正常则会减少再感染的概率

    Parameters
    ----------
    b : int
        抽样类型1的数量
    r : int
        抽样类型2的数量
    m : int
        抽到m个类型1
    n : int
        抽到n个类型2
    c : int
        每次抽完之后会增加c个同类型的数量

    Returns
    -------
    result : float
        算得的传染病模型的概率值

    Examples
    --------
    >>> infection(5, 6, 1, 2, 3)
    0.3093964858670741
    """

    result = comb(m+n,m)   # 从总数中哪几个是类型1选m个

    # 类型1的概率贡献
    for i in range(m):
        up = b + i*c
        down = b + r + i*c
        result *= up/down

    # 类型2的概率贡献
    for j in range(n):
        up = r + j*c
        # 注意这里的分母已经经过了m次类型1的累加
        down = b + r + (m+j)*c
        result *= up/down

    return result


def security(b:int, r:int, m:int, n:int, d:int, name1:str | int=1, name2:str | int=1, verbose=False) -> dict | float:
    """
    security(b, r, m, n, d, /) -> float \n
    每次抽样都会放入d个不同类型的样本，即检测安全，那么就放松，下一次的事故率就增加，
    检测事故，那么就加紧，下一次的事故率就下降  \n
    但是注意，这个模型没有通解，抽取的顺序影响概率的计算，因此需要先得到要计算的全排列
    Parameters
    ----------
    b : int
        抽样类型1的数量
    r : int
        抽样类型2的数量
    m : int
        抽到m个类型1
    n : int
        抽到n个类型2
    d : int
        每次抽完之后会增加d个不同类型的数量
    name1 : str | int
        类型1的名称
    name2 : str | int
        类型2的名称
    verbose : bool
        打印运算过程

    Returns
    -------
    out : dict | float
        如果verbose=True,计算不同组合下的安全模型的概率值
        如果verbose=False,直接输出所有组合下的概率和

    Examples
    --------
    >>> security(5, 6, 1, 2, 3, verbose=False)
    0.47097020626432384
    >>> security(5, 6, 1, 2, 3, verbose=True)
    {'1 0 1': 0.16501145912910617, '0 1 1': 0.15469824293353704, '1 1 0': 0.15126050420168066}


    每次抽完之后,放入3个不同类型的样本
    """
    import itertools
    name1, name2 = str(name1), str(name2)  #将名称统一转为字符串，便于后续比较和输出
    array = [name1] * m + [name2] * n
    result={}
    #生成所有不同排列
    # itertools.permutations————生成所有排列（含重复）  set()————去掉重复的排列
    permutations = set(itertools.permutations(array))

    for group in permutations:
        # 遍历所有排列
        _r = r
        _b = b
        out = 1
        for sample in group:
            if sample == name1:
                out *= _b / (_b + _r)
                _r += d
            elif sample == name2:
                out *= _r / (_b + _r)
                _b += d
        result[" ".join(group)] = out
    if verbose:
        return result
    else:
        return sum(result.vakues())
