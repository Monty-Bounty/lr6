% ЗАДАНИЕ:
% Уравнение: x*(2*x-1)*y' + y^2 - (4*x+1)*y + 4*x = 0
% y' = (-y^2 + (4*x+1)*y - 4*x) / (x*(2*x-1))
% Начальное условие: y(0) = 1
% Интервал: x от 0 до 1
% Шаг: h = 0.1

clear;
clc;
close all;

x0 = 0.0;
y0 = 1.0;
h = 0.1;
x_end = 1.0;

% Создаем вектор узлов сетки
x_points = x0:h:x_end;
n_points = length(x_points);

% Инициализация векторов для хранения решений
y_euler = zeros(1, n_points);
y_rk4 = zeros(1, n_points);

% Задаем начальные условия для каждого метода
y_euler(1) = y0;
y_rk4(1) = y0;

% Численное решение

% Метод Эйлера
for i = 1:(n_points - 1)
    x_current = x_points(i);
    y_current = y_euler(i);
    
    % Расчет производной y' = f(x,y)
    if abs(x_current - 0.0) < 1e-9
        dydx = -1.0;
    elseif abs(x_current - 0.5) < 1e-9
        dydx = 2.0 / 3.0;
    else
        num = -y_current^2 + (4*x_current + 1)*y_current - 4*x_current;
        dem = x_current * (2*x_current - 1);
        dydx = num / dem;
    end
    
    y_euler(i+1) = y_current + h * dydx;
end

% Метод Рунге-Кутта 4-го порядка
for i = 1:(n_points - 1)
    xi = x_points(i);
    yi = y_rk4(i);
    
    % Расчет k1
    x_k = xi;
    y_k = yi;
    if abs(x_k - 0.0) < 1e-9
        dydx1 = -1.0;
    elseif abs(x_k - 0.5) < 1e-9
        dydx1 = 2.0 / 3.0;
    else
        num = -y_k^2 + (4*x_k + 1)*y_k - 4*x_k;
        dem = x_k * (2*x_k - 1);
        dydx1 = num / dem;
    end
    k1 = h * dydx1;
    
    % Расчет k2
    x_k = xi + 0.5 * h;
    y_k = yi + 0.5 * k1;
    if abs(x_k - 0.0) < 1e-9
        dydx2 = -1.0;
    elseif abs(x_k - 0.5) < 1e-9
        dydx2 = 2.0 / 3.0;
    else
        num = -y_k^2 + (4*x_k + 1)*y_k - 4*x_k;
        dem = x_k * (2*x_k - 1);
        dydx2 = num / dem;
    end
    k2 = h * dydx2;

    % Расчет k3
    x_k = xi + 0.5 * h;
    y_k = yi + 0.5 * k2;
    if abs(x_k - 0.0) < 1e-9
        dydx3 = -1.0; 
    elseif abs(x_k - 0.5) < 1e-9
        dydx3 = 2.0 / 3.0;
    else
        num = -y_k^2 + (4*x_k + 1)*y_k - 4*x_k;
        dem = x_k * (2*x_k - 1);
        dydx3 = num / dem;
    end
    k3 = h * dydx3;
    
    % Расчет k4
    x_k = xi + h;
    y_k = yi + k3;
    if abs(x_k - 0.0) < 1e-9
        dydx4 = -1.0;
    elseif abs(x_k - 0.5) < 1e-9
        dydx4 = 2.0 / 3.0;
    else
        num = -y_k^2 + (4*x_k + 1)*y_k - 4*x_k;
        dem = x_k * (2*x_k - 1);
        dydx4 = num / dem;
    end
    k4 = h * dydx4;

    % Итоговый шаг
    y_rk4(i+1) = yi + (k1 + 2*k2 + 2*k3 + k4) / 6;
end

%Точное решение и вывод результатов
% функция для точного решения y(x) = (2*x^2 + 1) / (x + 1)
exact_solution = @(x) (2*x.^2 + 1) ./ (x + 1);

% Вычисляем точное решение в узлах сетки
y_exact = exact_solution(x_points);

fprintf('=========================================================================\n');
fprintf('|   x   |   Точное y(x)  |  Метод Эйлера  |  Метод Рунге-Кутта  |\n');
fprintf('=========================================================================\n');
for i = 1:n_points
    fprintf('| %5.1f |    %10.4f  |   %10.4f   |     %10.4f      |\n', ...
        x_points(i), y_exact(i), y_euler(i), y_rk4(i));
end
fprintf('=========================================================================\n');

%Построение графиков
figure('Name', 'Сравнение численных методов', 'NumberTitle', 'off'); 
hold on;

% Строим графики для каждого решения
plot(x_points, y_exact, 'o-', 'LineWidth', 2, 'DisplayName', 'Точное решение');
plot(x_points, y_euler, 's--', 'LineWidth', 1.5, 'DisplayName', 'Метод Эйлера');
plot(x_points, y_rk4, '^-', 'LineWidth', 1.5, 'DisplayName', 'Метод Рунге-Кутта 4');

% Настройка внешнего вида графика
title('Сравнение численных методов решения задачи Коши');
xlabel('x');
ylabel('y(x)');
legend('show', 'Location', 'northwest');
