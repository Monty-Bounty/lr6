clear;
clc;
close all;

% Определение функций системы ---
% Функции для правых частей системы дифференциальных уравнений
% y1' = y2
f1 = @(x, y1, y2) y2;
% y2' = -2*y1 - 3*y2 + 1 / (exp(x) + 1)
f2 = @(x, y1, y2) -2*y1 - 3*y2 + 1 / (exp(x) + 1);

x_start = 0;
x_end = 1;
h = 0.1;

% Начальные значения
y1_0 = 2 * (log(2) + 1);  % y(0)
y2_0 = -3 * log(2) - 1;   % y'(0)

% Создаем массив точек x, в которых будем искать решение
x_points = x_start:h:x_end;
num_steps = length(x_points);

% Инициализируем массивы для хранения решений
y1_rk4 = zeros(1, num_steps);
y2_rk4 = zeros(1, num_steps);

% Задаем начальные условия
y1_rk4(1) = y1_0;
y2_rk4(1) = y2_0;

%Реализация метода Рунге-Кутты 4-го порядка
for i = 1:(num_steps - 1)
    % Текущие значения x, y1, y2
    x_i = x_points(i);
    y1_i = y1_rk4(i);
    y2_i = y2_rk4(i);
    
    % Вычисление коэффициентов для системы
    k1 = h * f1(x_i, y1_i, y2_i);
    p1 = h * f2(x_i, y1_i, y2_i);
    
    k2 = h * f1(x_i + h/2, y1_i + k1/2, y2_i + p1/2);
    p2 = h * f2(x_i + h/2, y1_i + k1/2, y2_i + p1/2);
    
    k3 = h * f1(x_i + h/2, y1_i + k2/2, y2_i + p2/2);
    p3 = h * f2(x_i + h/2, y1_i + k2/2, y2_i + p2/2);
    
    k4 = h * f1(x_i + h, y1_i + k3, y2_i + p3);
    p4 = h * f2(x_i + h, y1_i + k3, y2_i + p3);
    
    % Расчет следующего значения y1 и y2
    y1_rk4(i+1) = y1_i + (k1 + 2*k2 + 2*k3 + k4) / 6;
    y2_rk4(i+1) = y2_i + (p1 + 2*p2 + 2*p3 + p4) / 6;
end

% Точное решение для сравнения
exact_solution = @(x) (exp(-x) + exp(-2*x)) .* (log(exp(x) + 1) + 1);
y_exact = exact_solution(x_points);

fprintf('Сравнение ручной реализации и точного решения в MATLAB \n');
fprintf('----------------------------------------------------------------\n');
fprintf('%-10s | %-25s | %-25s\n', 'x', 'y (Рунге-Кутта 4)', 'y (Точное)');
fprintf('----------------------------------------------------------------\n');
for i = 1:num_steps
    fprintf('%-10.4f | %-25.4f | %-25.4f\n', x_points(i), y1_rk4(i), y_exact(i));
end
fprintf('----------------------------------------------------------------\n');

% Построение графиков
figure('Name', 'Сравнение решений в MATLAB (ручная реализация)', 'NumberTitle', 'off');
hold on;

% График численного решения
plot(x_points, y1_rk4, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'Метод Рунге-Кутты 4');

% График точного решения
plot(x_points, y_exact, 'r--', 'LineWidth', 2, 'DisplayName', 'Точное решение');

% Настройка графика
title('Решение задачи Коши (ручная реализация)');
xlabel('x');
ylabel('y(x)');
legend('show', 'Location', 'southwest');
grid on;
hold off;

