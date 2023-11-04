import matplotlib.pyplot as plt


class SIR:
    population = 1_000_000
    immune = 0.0
    r0 = 1.2
    r0_lockdown = 0.5
    infected_duration = 2
    critical_duration = 20
    critical_ratio = 0.02

    lockdown_on_threshold = 0.005
    lockdown_off_threshold = 0.001

    def evolution(self, initial_infected):
        susceptible = self.population - initial_infected - self.immune
        infected = initial_infected
        resolved = self.immune
        critical = 0.0

        beta = self.r0 / self.infected_duration
        gamma = 1.0 / self.infected_duration
        delta = 1.0 / self.critical_duration

        sequence = [(susceptible, infected, resolved, critical)]

        while True:
            if infected / self.population >= self.lockdown_on_threshold:
                beta = self.r0_lockdown / self.infected_duration
            if infected / self.population < self.lockdown_off_threshold:
                beta = self.r0 / self.infected_duration

            newly_infected = beta * susceptible / self.population * infected
            newly_resolved = gamma * infected

            susceptible -= newly_infected
            infected += newly_infected - newly_resolved
            resolved += newly_resolved
            critical += self.critical_ratio * newly_resolved - delta * critical

            sequence += [(susceptible, infected, resolved, critical)]

            if infected < 1:
                break

        return sequence

    def plot(self, sequence):
        population = [self.population] * len(sequence)
        susceptible = [s for (s, i, r, c) in sequence]
        infected = [i for (s, i, r, c) in sequence]
        resolved = [r for (s, i, r, c) in sequence]
        critical = [c for (s, i, r, c) in sequence]

        plt.grid()
        # plt.plot(susceptible)
        plt.plot(infected)
        # plt.plot(resolved)
        plt.plot(critical)
        # plt.plot(population)
        plt.show()


m = SIR()
m.population = 2_090_000

seq = m.evolution(100)

print(
    "Attack: {:.3g}% ({:.3g}%), Critical {:.0f} (max {:.0f})".format(
        (seq[-1][2] - m.immune) / m.population * 100,
        seq[-1][2] / m.population * 100,
        (seq[-1][2] - m.immune) * m.critical_ratio,
        max(c for (s, i, r, c) in seq),
    )
)
m.plot(seq)
