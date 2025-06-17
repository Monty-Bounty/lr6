import numpy as np
import matplotlib.pyplot as plt

# --- 1. Определение системы дифференциальных уравнений ---
# Исходное уравнение: y'' + 3y' + 2y = 1 / (e^x + 1)
# Производим замену:
# y1 = y
# y2 = y'
# Получаем систему уравнений первого порядка:
# y1' = y2
# y2' = y'' = -2y1 - 3y2 + 1 / (e^x + 1)

def f1(x, y1, y2):
    """ y1' = y2 """
    return y2

def f2(x, y1, y2):
    """ y2' = -2*y1 - 3*y2 + 1 / (exp(x) + 1) """
    return -2 * y1 - 3 * y2 + 1 / (np.exp(x) + 1)

# --- 2. Точное решение ---
def exact_solution(x):
    """ y(x) = (e^(-x) + e^(-2*x)) * [ln(e^x + 1) + 1] """
    return (np.exp(-x) + np.exp(-2*x)) * (np.log(np.exp(x) + 1) + 1)

# --- 3. Начальные условия и параметры ---
x0 = 0.0
y1_0 = 2 * (np.log(2) + 1)  # y(0)
y2_0 = -3 * np.log(2) - 1   # y'(0)

# Интервал и шаг
x_end = 1.0
h = 0.1

# Массив точек x, в которых будем искать решение
x_points = np.arange(x0, x_end + h, h)
n_steps = len(x_points)

# --- 4. Метод Эйлера ---
y1_euler = np.zeros(n_steps)
y2_euler = np.zeros(n_steps)
y1_euler[0] = y1_0
y2_euler[0] = y2_0

for i in range(n_steps - 1):
    x_i = x_points[i]
    y1_i = y1_euler[i]
    y2_i = y2_euler[i]

    # Расчет по формулам Эйлера для системы
    y1_euler[i+1] = y1_i + h * f1(x_i, y1_i, y2_i)
    y2_euler[i+1] = y2_i + h * f2(x_i, y1_i, y2_i)

# --- 5. Метод Рунге-Кутты 4-го порядка ---
y1_rk4 = np.zeros(n_steps)
y2_rk4 = np.zeros(n_steps)
y1_rk4[0] = y1_0
y2_rk4[0] = y2_0

for i in range(n_steps - 1):
    x_i = x_points[i]
    y1_i = y1_rk4[i]
    y2_i = y2_rk4[i]

    # Расчет коэффициентов для y1
    k1 = h * f1(x_i, y1_i, y2_i)
    p1 = h * f2(x_i, y1_i, y2_i)

    k2 = h * f1(x_i + h/2, y1_i + k1/2, y2_i + p1/2)
    p2 = h * f2(x_i + h/2, y1_i + k1/2, y2_i + p1/2)

    k3 = h * f1(x_i + h/2, y1_i + k2/2, y2_i + p2/2)
    p3 = h * f2(x_i + h/2, y1_i + k2/2, y2_i + p2/2)

    k4 = h * f1(x_i + h, y1_i + k3, y2_i + p3)
    p4 = h * f2(x_i + h, y1_i + k3, y2_i + p3)

    # Обновление значений y1 и y2
    y1_rk4[i+1] = y1_i + (k1 + 2*k2 + 2*k3 + k4) / 6
    y2_rk4[i+1] = y2_i + (p1 + 2*p2 + 2*p3 + p4) / 6

# --- 6. Вывод результатов ---

# Вывод в консоль
print("--- Сравнение численных методов с точным решением ---")
print("-" * 70)
print(f"{'x':<10} | {'y (Эйлер)':<15} | {'y (Рунге-Кутта)':<20} | {'y (Точное)':<15}")
print("-" * 70)

y_exact_points = exact_solution(x_points)

for i in range(n_steps):
    print(f"{x_points[i]:<10.2f} | {y1_euler[i]:<15.4f} | {y1_rk4[i]:<20.4f} | {y_exact_points[i]:<15.4f}")

print("-" * 70)

# Построение графиков
plt.style.use('seaborn-v0_8-whitegrid')
plt.figure(figsize=(12, 8))
plt.plot(x_points, y_exact_points, 'k-', label='Точное решение', linewidth=3)
plt.plot(x_points, y1_euler, 'r--o', label='Метод Эйлера', markersize=5)
plt.plot(x_points, y1_rk4, 'b--^', label='Метод Рунге-Кутты 4-го порядка', markersize=5)

plt.title('Решение задачи Коши', fontsize=16)
plt.xlabel('x', fontsize=12)
plt.ylabel('y(x)', fontsize=12)
plt.legend(fontsize=12)
plt.grid(True)
plt.show()


with open("results.txt", "w", encoding="utf-8") as f:
    f.write("--- Сравнение численных методов с точным решением ---\n")
    f.write("-" * 70 + "\n")
    f.write(f"{'x':<10} | {'y (Эйлер)':<15} | {'y (Рунге-Кутта)':<20} | {'y (Точное)':<15}\n")
    f.write("-" * 70 + "\n")
    for i in range(n_steps):
        f.write(f"{x_points[i]:<10.2f} | {y1_euler[i]:<15.4f} | {y1_rk4[i]:<20.4f} | {y_exact_points[i]:<15.4f}\n")
    f.write("-" * 70 + "\n")
print("\nРезультаты также сохранены в файл 'results.txt'")