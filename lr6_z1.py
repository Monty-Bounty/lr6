import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# ЗАДАНИЕ:
# Уравнение: x*(2*x-1)*y' + y^2 - (4*x+1)*y + 4*x = 0
# Начальное условие: y(0) = 1
# Интервал: x от 0 до 1
# Шаг: h = 0.1
#
# Для использования численных методов, выразим y' (производную):
# y' = (-y^2 + (4*x+1)*y - 4*x) / (x*(2*x-1))
# ---------------------------------------------------------------------------

def f(x, y):
    """
    Правая часть дифференциального уравнения y' = f(x, y).
    """
    # ОСОБЕННОСТЬ 1: Устранимый разрыв в x=0.
    # В точке x=0 знаменатель обращается в ноль. Предел y' при x -> 0 равен -1.
    if np.isclose(x, 0.0):
        return -1.0
    
    # ОСОБЕННОСТЬ 2: Устранимый разрыв в x=0.5.
    # В точке x=0.5 знаменатель также обращается в ноль.
    # Можно показать (например, с помощью правила Лопиталя), что предел
    # y' при x -> 0.5 для точного решения равен 2/3.
    if np.isclose(x, 0.5):
        return 2.0 / 3.0

    denominator = x * (2 * x - 1)
    numerator = -y**2 + (4 * x + 1) * y - 4 * x
    return numerator / denominator

def exact_solution(x):
    """
    Точное решение задачи Коши для сравнения.
    y(x) = (2*x^2 + 1) / (x + 1)
    Эта функция теперь векторизована и работает с массивами NumPy.
    """
    return (2 * np.power(x, 2) + 1) / (x + 1)

# --- Параметры задачи ---
x0 = 0.0
y0 = 1.0
h = 0.1
x_end = 1.0

# Создаем массив точек x, на которых будем искать решение
x_points = np.arange(x0, x_end + h, h)
n_points = len(x_points)

# --- Инициализация массивов для хранения решений ---
y_euler = np.zeros(n_points)
y_rk4 = np.zeros(n_points)
y_adams = np.zeros(n_points)

# Задаем начальные условия
y_euler[0] = y0
y_rk4[0] = y0
y_adams[0] = y0

# --- 1. Метод Эйлера ---
for i in range(n_points - 1):
    y_euler[i+1] = y_euler[i] + h * f(x_points[i], y_euler[i])

# --- 2. Метод Рунге-Кутта 4-го порядка ---
for i in range(n_points - 1):
    xi = x_points[i]
    yi = y_rk4[i]
    
    k1 = h * f(xi, yi)
    k2 = h * f(xi + 0.5 * h, yi + 0.5 * k1)
    k3 = h * f(xi + 0.5 * h, yi + 0.5 * k2)
    k4 = h * f(xi + h, yi + k3)
    
    y_rk4[i+1] = yi + (k1 + 2*k2 + 2*k3 + k4) / 6


y_exact = exact_solution(x_points)

# --- Вывод результатов в консоль в виде таблицы ---
print("="*65)
print(f"{'x':>5} | {'Точное y(x)':>15} | {'Метод Эйлера':>15} | {'Метод Рунге-Кутта':>20}")
print("-"*65)
for i in range(n_points):
    # Используем np.nan_to_num, чтобы заменить возможные nan/inf на 0 для вывода,
    # хотя после исправлений их быть не должно.
    print(f"{x_points[i]:>5.1f} | {y_exact[i]:>15.4f} | {np.nan_to_num(y_euler[i]):>15.4f} | {np.nan_to_num(y_rk4[i]):>20.4f} ")
print("="*65)


# --- Построение графиков ---
plt.style.use('seaborn-v0_8-whitegrid')
plt.figure(figsize=(12, 8))

plt.plot(x_points, y_exact, 'o-', label='Точное решение', linewidth=2, zorder=5)
plt.plot(x_points, y_euler, 's--', label='Метод Эйлера', alpha=0.8)
plt.plot(x_points, y_rk4, '^-', label='Метод Рунге-Кутта 4', alpha=0.8)

# Настройка графика
plt.title('Сравнение численных методов решения задачи Коши', fontsize=16)
plt.xlabel('x', fontsize=12)
plt.ylabel('y(x)', fontsize=12)
plt.legend(fontsize=12)
plt.grid(True)
plt.show()
