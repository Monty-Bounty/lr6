import numpy as np
import matplotlib.pyplot as plt

# Система дифференциальных уравнений:
# y1' = 5*y1 - 3*y2 + 2*e^(3x)
# y2' = y1 + y2 + 5*e^(-x)

def f1(x, y1, y2):
  return 5 * y1 - 3 * y2 + 2 * np.exp(3 * x)

def f2(x, y1, y2):
  return y1 + y2 + 5 * np.exp(-x)

y1_0 = -1.0
y2_0 = -2.0

x_start = 0.0
x_end = 1.0
h = 0.1

def exact_y1(x):
  return np.exp(2 * x) + 3 * np.exp(4 * x) - np.exp(-x) - 4 * np.exp(3 * x)

def exact_y2(x):
  return np.exp(2 * x) + np.exp(4 * x) - 2 * np.exp(-x) - 2 * np.exp(3 * x)

def euler_method():
  """
  Решает систему ОДУ методом Эйлера.
  Возвращает списки со значениями x, y1 и y2.
  """
  # Создаем массивы для хранения результатов
  n_steps = int((x_end - x_start) / h) + 1
  x_values = np.linspace(x_start, x_end, n_steps)
  y1_values = np.zeros(n_steps)
  y2_values = np.zeros(n_steps)

  # Задаем начальные условия
  y1_values[0] = y1_0
  y2_values[0] = y2_0

  # Итеративно вычисляем решение
  for i in range(n_steps - 1):
    x_i = x_values[i]
    y1_i = y1_values[i]
    y2_i = y2_values[i]

    y1_values[i + 1] = y1_i + h * f1(x_i, y1_i, y2_i)
    y2_values[i + 1] = y2_i + h * f2(x_i, y1_i, y2_i)

  return x_values, y1_values, y2_values

def runge_kutta_4_method():
  """
  Решает систему ОДУ методом Рунге-Кутта 4-го порядка.
  Возвращает списки со значениями x, y1 и y2.
  """
  # Создаем массивы для хранения результатов
  n_steps = int((x_end - x_start) / h) + 1
  x_values = np.linspace(x_start, x_end, n_steps)
  y1_values = np.zeros(n_steps)
  y2_values = np.zeros(n_steps)

  # Задаем начальные условия
  y1_values[0] = y1_0
  y2_values[0] = y2_0

  # Итеративно вычисляем решение
  for i in range(n_steps - 1):
    x_i = x_values[i]
    y1_i = y1_values[i]
    y2_i = y2_values[i]

    # Вычисляем коэффициенты для y1
    k1 = h * f1(x_i, y1_i, y2_i)
    p1 = h * f2(x_i, y1_i, y2_i)

    k2 = h * f1(x_i + h / 2, y1_i + k1 / 2, y2_i + p1 / 2)
    p2 = h * f2(x_i + h / 2, y1_i + k1 / 2, y2_i + p1 / 2)

    k3 = h * f1(x_i + h / 2, y1_i + k2 / 2, y2_i + p2 / 2)
    p3 = h * f2(x_i + h / 2, y1_i + k2 / 2, y2_i + p2 / 2)

    k4 = h * f1(x_i + h, y1_i + k3, y2_i + p3)
    p4 = h * f2(x_i + h, y1_i + k3, y2_i + p3)

    # Вычисляем следующие значения y1 и y2
    y1_values[i + 1] = y1_i + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    y2_values[i + 1] = y2_i + (p1 + 2 * p2 + 2 * p3 + p4) / 6

  return x_values, y1_values, y2_values


# Выполняем вычисления
x_euler, y1_euler, y2_euler = euler_method()
x_rk4, y1_rk4, y2_rk4 = runge_kutta_4_method()

# Вычисляем точные значения для тех же точек x
y1_exact_vals = exact_y1(x_rk4)
y2_exact_vals = exact_y2(x_rk4)

print("--- Результаты численного решения ---")
print("-" * 110)
print(f"{'x':>5} | {'y1 (Эйлер)':>15} | {'y2 (Эйлер)':>15} | {'y1 (Р-К 4)':>15} | {'y2 (Р-К 4)':>15} | {'y1 (Точное)':>15} | {'y2 (Точное)':>15}")
print("-" * 110)

for i in range(len(x_rk4)):
  print(f"{x_rk4[i]:>5.2f} | {y1_euler[i]:>15.4f} | {y2_euler[i]:>15.4f} | {y1_rk4[i]:>15.4f} | {y2_rk4[i]:>15.4f} | {y1_exact_vals[i]:>15.4f} | {y2_exact_vals[i]:>15.4f}")

print("-" * 110)

with open('results.txt', 'w', encoding='utf-8') as f:
  f.write("--- Результаты численного решения ---\n")
  f.write("-" * 110 + "\n")
  f.write(f"{'x':>5} | {'y1 (Эйлер)':>15} | {'y2 (Эйлер)':>15} | {'y1 (Р-К 4)':>15} | {'y2 (Р-К 4)':>15} | {'y1 (Точное)':>15} | {'y2 (Точное)':>15}\n")
  f.write("-" * 110 + "\n")

  for i in range(len(x_rk4)):
    f.write(f"{x_rk4[i]:>5.2f} | {y1_euler[i]:>15.4f} | {y2_euler[i]:>15.4f} | {y1_rk4[i]:>15.4f} | {y2_rk4[i]:>15.4f} | {y1_exact_vals[i]:>15.4f} | {y2_exact_vals[i]:>15.4f}\n")
  
  f.write("-" * 110 + "\n")

print("\nРезультаты также сохранены в файл 'results.txt'")

plt.style.use('seaborn-v0_8-whitegrid')
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
fig.suptitle('Сравнение численных методов с точным решением', fontsize=16)

# График для y1(x)
ax1.plot(x_euler, y1_euler, 'o--', label='Метод Эйлера', color='orange')
ax1.plot(x_rk4, y1_rk4, 's-', label='Метод Рунге-Кутта 4', color='blue')
ax1.plot(x_rk4, y1_exact_vals, '-', label='Точное решение', color='black', linewidth=2)
ax1.set_ylabel('y1(x)', fontsize=12)
ax1.legend()
ax1.set_title('Решение для y1(x)', fontsize=14)

# График для y2(x)
ax2.plot(x_euler, y2_euler, 'o--', label='Метод Эйлера', color='orange')
ax2.plot(x_rk4, y2_rk4, 's-', label='Метод Рунге-Кутта 4', color='blue')
ax2.plot(x_rk4, y2_exact_vals, '-', label='Точное решение', color='black', linewidth=2)
ax2.set_xlabel('x', fontsize=12)
ax2.set_ylabel('y2(x)', fontsize=12)
ax2.legend()
ax2.set_title('Решение для y2(x)', fontsize=14)

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()



'''
    x | y1 (Эйлер) | y2 (Эйлер) |    y1 (Р-К 4) |    y2 (Р-К 4) | y1 (Точное) | y2 (Точное)
'''