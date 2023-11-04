#!/usr/bin/python

import math
import matplotlib.pyplot as plt


def ramanujan_comb(n):
    return (
        n * math.log(n)
        - n
        + (math.log(n * (1 + 4 * n * (1 + 2 * n)))) / 6
        + math.log(math.pi) / 2
    )


# Probability for at least i of of n events ocurring
def p_le(n, i, p):
    return sum(math.comb(n, j) * p ** j * (1 - p) ** (n - j) for j in range(i, n + 1))

print(p_le(600, 300, 0.15))
