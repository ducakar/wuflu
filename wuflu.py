import math
import datetime as dt
import matplotlib.pyplot as plt
import json


def lin_fit(yy):
    """Aproksimacija z linearno funkcijo po metodi najmanjših kvadratov."""
    """https://mathworld.wolfram.com/LeastSquaresFitting.html"""
    """y = a + b * x"""
    n = len(yy)
    xx = [x for x in range(n)]

    sum_x = sum(xx)
    sum_y = sum(yy)
    sum_x_2 = sum(x ** 2 for x in xx)
    sum_x_y = sum(x * y for (x, y) in zip(xx, yy))

    a = (sum_y * sum_x_2 - sum_x * sum_x_y) / (n * sum_x_2 - sum_x ** 2)
    b = (n * sum_x_y - sum_x * sum_y) / (n * sum_x_2 - sum_x ** 2)
    return a, b


def exp_fit(yy):
    """Aproksimacija z eksponentno funkcijo po metodi najmanjših kvadratov."""
    """https: // mathworld.wolfram.com/LeastSquaresFittingExponential.html"""
    """y = a * b^x"""
    n = len(yy)
    xx = [x for x in range(n)]

    sum_x = sum(xx)
    sum_x_2 = sum(x ** 2 for x in xx)
    sum_log_y = sum(math.log(y) for y in yy)
    sum_x_log_y = sum(x * math.log(y) for (x, y) in zip(xx, yy))

    a = (sum_log_y * sum_x_2 - sum_x * sum_x_log_y) / (n * sum_x_2 - sum_x ** 2)
    b = (n * sum_x_log_y - sum_x * sum_log_y) / (n * sum_x_2 - sum_x ** 2)
    return math.exp(a), math.exp(b)


with open("data.json") as f:
    data = json.load(f)

infections = [x["cases"]["active"] for x in data]

days = [dt.date(2020, 3, 4) + dt.timedelta(t) for t, _ in enumerate(infections)]
infections = infections[10:25]
days = days[10:25]

lin_a, lin_b = lin_fit(infections)
exp_a, exp_b = exp_fit(infections)
predicted_lin = [lin_a + lin_b * t for t, _ in enumerate(infections)]
predicted_exp = [exp_a * exp_b ** t for t, _ in enumerate(infections)]

print("{} \\times {}^n".format(exp_a, exp_b))

_, ax = plt.subplots()
ax.set_xticks(days)
ax.set_xticklabels(["{}.{}.".format(d.day, d.month) for d in days])
# ax.set_yscale("log")
plt.plot(days, infections)
plt.plot(days, predicted_lin)
plt.plot(days, predicted_exp)
plt.scatter(days, infections)
plt.grid(axis="y")
plt.show()
