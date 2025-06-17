clear;
clc;
close all;

% y1' = 5*y1 - 3*y2 + 2*exp(3x)
% y2' = y1 + y2 + 5*exp(-x)

f1 = @(x, y1, y2) 5*y1 - 3*y2 + 2*exp(3*x);
f2 = @(x, y1, y2) y1 + y2 + 5*exp(-x);

% Точное решение для сравнения
y1_exact_fun = @(x) exp(2*x) + 3*exp(4*x) - exp(-x) - 4*exp(3*x);
y2_exact_fun = @(x) exp(2*x) + exp(4*x) - 2*exp(-x) - 2*exp(3*x);


x_start = 0.0;
x_end = 1.0;
h = 0.1;
y0 = [-1; -2];

% Сетка для вычислений
x_span = x_start:h:x_end;
n_steps = length(x_span);

% Метод Эйлера
y_euler = zeros(2, n_steps); % Массив для результатов [y1; y2]
y_euler(:, 1) = y0;          % Записываем начальные условия

for i = 1:(n_steps - 1)
    x_i = x_span(i);
    y1_i = y_euler(1, i);
    y2_i = y_euler(2, i);
    
    % Вычисляем следующие значения по формулам Эйлера
    y_euler(1, i+1) = y1_i + h * f1(x_i, y1_i, y2_i);
    y_euler(2, i+1) = y2_i + h * f2(x_i, y1_i, y2_i);
end

% Метод Рунге-Кутта 4-го порядка
y_rk4 = zeros(2, n_steps);   % Массив для результатов [y1; y2]
y_rk4(:, 1) = y0;            % Записываем начальные условия

for i = 1:(n_steps - 1)
    x_i = x_span(i);
    y1_i = y_rk4(1, i);
    y2_i = y_rk4(2, i);
    
    % Вычисляем коэффициенты для y1 и y2
    k1 = h * f1(x_i, y1_i, y2_i);
    p1 = h * f2(x_i, y1_i, y2_i);
    
    k2 = h * f1(x_i + h/2, y1_i + k1/2, y2_i + p1/2);
    p2 = h * f2(x_i + h/2, y1_i + k1/2, y2_i + p1/2);
    
    k3 = h * f1(x_i + h/2, y1_i + k2/2, y2_i + p2/2);
    p3 = h * f2(x_i + h/2, y1_i + k2/2, y2_i + p2/2);
    
    k4 = h * f1(x_i + h, y1_i + k3, y2_i + p3);
    p4 = h * f2(x_i + h, y1_i + k3, y2_i + p3);
    
    % Вычисляем следующие значения
    y_rk4(1, i+1) = y1_i + (k1 + 2*k2 + 2*k3 + k4) / 6;
    y_rk4(2, i+1) = y2_i + (p1 + 2*p2 + 2*p3 + p4) / 6;
end

% Вычисление точных значений
y1_exact_vals = y1_exact_fun(x_span);
y2_exact_vals = y2_exact_fun(x_span);

% Вывод результатов в командное окно
fprintf('Результаты численного решения\n');
fprintf('--------------------------------------------------------------------------------------------\n');
fprintf('%5s | %15s | %15s | %15s | %15s | %15s | %15s\n', 'x', 'y1 (Эйлер)', 'y2 (Эйлер)', 'y1 (Р-К 4)', 'y2 (Р-К 4)', 'y1 (Точное)', 'y2 (Точное)');
fprintf('--------------------------------------------------------------------------------------------\n');

for i = 1:n_steps
    fprintf('%5.2f | %10.4f | %10.4f | %10.4f | %10.4f | %10.4f | %10.4f\n', ...
        x_span(i), y_euler(1,i), y_euler(2,i), y_rk4(1,i), y_rk4(2,i), y1_exact_vals(i), y2_exact_vals(i));
end
fprintf('--------------------------------------------------------------------------------------------\n');


% Построение графиков
figure('Name', 'Сравнение численных методов', 'NumberTitle', 'off');

% График для y1(x)
subplot(2, 1, 1);
hold on;
plot(x_span, y_euler(1,:), 'o--', 'DisplayName', 'Метод Эйлера', 'Color', [0.9290 0.6940 0.1250]);
plot(x_span, y_rk4(1,:), 's-', 'DisplayName', 'Метод Рунге-Кутта 4', 'Color', [0 0.4470 0.7410]);
plot(x_span, y1_exact_vals, '-', 'DisplayName', 'Точное решение', 'Color', 'black', 'LineWidth', 2);
title('Решение для y1(x)');
ylabel('y1(x)');
legend('show', 'Location', 'northwest');
grid on;
box on;
hold off;

% График для y2(x)
subplot(2, 1, 2);
hold on;
plot(x_span, y_euler(2,:), 'o--', 'DisplayName', 'Метод Эйлера', 'Color', [0.9290 0.6940 0.1250]);
plot(x_span, y_rk4(2,:), 's-', 'DisplayName', 'Метод Рунге-Кутта 4', 'Color', [0 0.4470 0.7410]);
plot(x_span, y2_exact_vals, '-', 'DisplayName', 'Точное решение', 'Color', 'black', 'LineWidth', 2);
title('Решение для y2(x)');
xlabel('x');
ylabel('y2(x)');
legend('show', 'Location', 'northwest');
grid on;
box on;
hold off;