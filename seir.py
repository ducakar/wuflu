# %%
import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = [6, 3.5]
plt.rcParams["figure.dpi"] = 96

# %%
class Lockdown:
    def __init__(self, R0, on_threshold, off_threshold):
        self.R0 = R0
        self.on_threshold = on_threshold
        self.off_threshold = off_threshold


class SEIR:
    P = 1.0
    R0 = 2.28
    chi = 0.05
    kappa = 0.01
    mu = 0.0066
    t_E = 4.5
    t_I = 2.0
    t_H = 25.0
    t_C = 28.0

    def __init__(self, title):
        self.title = title

    def evolution(self, infected, immune, lockdown=None):
        S = [self.P - infected - immune]
        E = [0.0]
        I = [infected]
        R = [immune]
        H = [0.0]
        C = [0.0]
        F = [0.0]
        R0 = self.R0

        while True:
            if lockdown:
                if I[-1] / self.P >= lockdown.on_threshold:
                    R0 = lockdown.R0
                elif I[-1] / self.P < lockdown.off_threshold:
                    R0 = self.R0

            dEp = R0 * S[-1] / self.P * I[-1] / self.t_I
            dIp = E[-1] / self.t_E
            dRp = I[-1] / self.t_I
            dS = -dEp
            dE = dEp - dIp
            dI = dIp - dRp
            dR = dRp
            dH = self.chi * dRp - H[-1] / self.t_H
            dC = self.kappa * dRp - C[-1] / self.t_C
            dF = self.mu * dRp

            S += [S[-1] + dS]
            E += [E[-1] + dE]
            I += [I[-1] + dI]
            R += [R[-1] + dR]
            H += [H[-1] + dH]
            C += [C[-1] + dC]
            F += [F[-1] + dF]

            if max(I[-1], H[-1], C[-1]) / self.P < 1e-6:
                break

        return S, E, I, R, H, C, F

    def plot(self, S, E, I, R, H, C, F):
        attack = (S[0] - S[-1]) / self.P
        r = 1.0 - S[-1] / self.P
        all_H = (S[0] - S[-1]) * self.chi
        max_H = max(H)
        all_C = (S[0] - S[-1]) * self.kappa
        max_C = max(C)
        all_F = F[-1]

        plt.title(self.title)
        plt.ylabel("populacija")
        plt.xlabel("dnevi")
        plt.grid()
        plt.plot(E, label="inkubacija")
        plt.plot(I, label="kužni")
        plt.plot(H, label="resni")
        # plt.plot(C, label="kritični", color=[1.0, 0.0, 0.0])
        plt.plot(F, label="mrtvi", color=[0, 0, 0])
        plt.legend()
        plt.show()

        print(
            "prekuženi %.3g%% (novi %.3g%%), resni %0.f (maks %.0f), kritični %.0f (maks %.0f), mrtvi %.0f"
            % (r * 100, attack * 100, all_H, max_H, all_C, max_C, all_F)
        )


# %%
m = SEIR("Covid-19")
m.P = 2_090_000
m.R0 = 2.28
S, E, I, R, H, C, F = m.evolution(100.0, 0.0)
m.plot(S, E, I, R, H, C, F)
