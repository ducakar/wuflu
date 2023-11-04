import json

with open("deaths-si.json") as f:
    monthly_deaths = json.load(f)
    monthly_average = []

for i in range(0, 12):
    s = sum(monthly_deaths[i + j * 12] for j in range(15, 20))
    monthly_average.append(round(s / 5))

for i, v in enumerate(monthly_deaths[20 * 12 :]):
    print(1 + i, v / monthly_average[i])
