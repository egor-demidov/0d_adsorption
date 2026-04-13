import numpy as np
import matplotlib.pyplot as plt


tau_1 = 1.0
tau_2 = 8.5
t_ads_start = 50.0
t_ads_end = 200.0

k_1 = 1.0/tau_1
k_2 = 1.0/tau_2
a = k_2 / (k_1 + k_2)
b = k_1 / (k_1 + k_2)
s = 1.0/k_2 * np.arctanh((k_1 - k_2)/(2.0 * k_1))

def g_of_t(t):
    if t < -s:
        return a + a * np.tanh(k_1 * (t + s))
    return a + b * np.tanh(k_2 * (t + s))

def f_of_t(t):
    return g_of_t(t - t_ads_start) - g_of_t(t - t_ads_end)

t_span = np.linspace(0.0, t_ads_end + 100.0, 500)
y_span = np.zeros(t_span.shape)

for i in range(len(t_span)):
    y_span[i] = f_of_t(t_span[i])

plt.plot(t_span, y_span)
plt.show()
