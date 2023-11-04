# %%
import math
import matplotlib.pyplot as plt

# %%
def normal(mu, sigma, x):
    return (
        1
        / (sigma * math.sqrt(2 * math.pi))
        * math.exp(-1 / 2 * (x - mu) ** 2 / sigma ** 2)
    )


# %%
xx = [i / 1000 for i in range(-3000, +3001)]
s1 = [normal(0.0, 0.5, x) for x in xx]
s2 = [normal(-0.05, 0.6, x) for x in xx]
s3 = [b / (a + b) for (a, b) in zip(s1, s2)]
plt.plot(s1)
plt.plot(s2)
plt.plot(s3)
plt.plot([0.5 for x in range(len(s3))])
plt.show()
