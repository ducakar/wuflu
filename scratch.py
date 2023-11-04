# %%
P = 2_100_000
I = 1_000
N = 10000

p = 1.0
for i in range(N):
    p *= (P - I - i) / (P - i)
    if i % 100 == 99:
        print("%d: %.2g%%" % (i + 1, (1 - p) * 100))
