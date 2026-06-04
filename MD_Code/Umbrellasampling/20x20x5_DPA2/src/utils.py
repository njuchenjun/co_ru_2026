import random
import math

def gaussian_knuth():
    """原 Fortran 中的高斯随机数生成器（Knuth 算法）"""
    a1 = 3.949846138
    a3 = 0.252408784
    a5 = 0.076542912
    a7 = 0.008355968
    a9 = 0.029899776
    uni = [random.random() for _ in range(12)]
    rrr = (sum(uni) - 6.0) / 4.0
    rr2 = rrr * rrr
    return rrr * (a1 + rr2 * (a3 + rr2 * (a5 + rr2 * (a7 + rr2 * a9))))

def fast_gaussian():
    """使用 Python 内置的高斯分布（更快）"""
    return random.gauss(0.0, 1.0)