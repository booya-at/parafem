import matplotlib.pyplot as plt
import numpy as np


def damping(t0, x0, f, h):
    return t0 + h, x0 - f * h * x0


def test(tend, f=1, h=0.01):
    t = [0]
    x = [1]
    for i in range(int(tend / h)):
        ti, xi = damping(t[-1], x[-1], f, h)
        t.append(ti)
        x.append(xi)
    return t, x

f = 1
tend = 10

plt.plot(*test(tend, f, h=0.01))
t = np.linspace(0, tend, 100)
plt.plot(t, np.e ** (-f * t))
plt.show()

