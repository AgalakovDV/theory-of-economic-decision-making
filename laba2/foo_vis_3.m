function [] = foo_vis_3(n)

format compact;
p = 3;
mu = zeros(1,p);
% mu = ones(1,p);
% for j = 1:p
%     mu(j) = j^5*(-1)^j;
% end;
Sigma = [ 10 5 30;
          5 50 10;
          30 10 100; ];
% Sigma = [ 451.39 271.17 168.70
%           271.17 171.73 103.29
%           168.70 103.29 66.65  ];
set = mvnrnd(mu, Sigma, n);
%центрируем данные
% mean_vec = mean(set);
% for i = 1:p
%     set(:,i) = set(:,i) - mean_vec(i);
% end;

% Приступим к переходу от x к z
Z = zeros(n, p);
% z_i = l_i*x, l_i - с.в. матрицы ковариаций вектора x

% Найдём собственные значения матрицы ковариаций
[L, R] = eig(Sigma);

% Отсортируем с.ч.
% r_sort - массив, отсортированных с.ч.
% ind - массив номеров, в каком порядке они находились.
[r_sort, ind] = sort(diag(R));
% disp('r_sort');
% disp(r_sort);
for k = 1:n
    for i = 1:p
        Z(k, i) = set(k,:) * L(:, ind(p + 1 - i));%возможно транспонировать
    end;
end;

% Найдем матрицу ковариаций через функцию cov
 Sigma_z1  = cov(Z); 

 % Найдём матрицу ковариаций по z
 Sigma_z2 = eye(p);
 for i = 1:p
     Sigma_z2(i,i) = r_sort(p + 1 - i);
 end;
% disp('set:');
% disp(set);
% disp('Z:');
% disp(Z);
% disp('L');
% disp(L);
% disp('R');
% disp(R);

% % % % % % % % % % %  disp('Sigma');
% % % % % % % % % % %  disp(Sigma);
% % % % % % % % % % %  disp('Sigma2');
% % % % % % % % % % %  disp(cov(set));
% disp('disp');
% disp(cov(set(:,1)));
% % % % % % % % % % %  disp('Sigma_z1')
% % % % % % % % % % %  disp(Sigma_z1);
disp('Sigma_z2')
disp(Sigma_z2);
 
% Изобразим точки X
figure('Name','Measured Data');
subplot(3,1,1);
title('Set');
hold on;
plot3(set(:,1), set(:,2), set(:,3), 'k.');
xlabel('x');
ylabel('y');
zlabel('z');
grid on;
view(50, 40);
subplot(3,1,2);
title('Z');
hold on;
plot3(Z(:,1), Z(:,2), Z(:,3), 'k.');
xlabel('x');
ylabel('y');
zlabel('z');
grid on;
view(50, 40);
subplot(3,1,3);
hold on;
title('Z_2');



% sum_0 = -1;
% % % % sum_disp = sum(diag(Sigma));
% disp('sum_disp');
% disp(sum_disp);
% Процентная доля дисперсии на каждую из компонент

% Sigma_proc - массив долей дисперсий для каждого измерения из p штук

Sigma_proc = zeros(1,p);
sum0 = 0;
for i = 1:p
    Sigma_proc(i) = sum0 + Sigma_z2(i,i) / sum(r_sort);
    sum0 = Sigma_proc(i);
end;

disp('Sigma_proc');
disp(Sigma_proc);

% Порог 70%
isImportant = zeros(1,p);
for i = 1:p
    isImportant(i) = 1;
    if (Sigma_proc(i) >= 0.7)
        break;
    end;
end;

% Создадим матрицу Z' из Z путем конкатенации столбцов 

subplot(3,1,3);
hold on;
if (sum(isImportant) == 2)
    if (isImportant(2) == isImportant(3))
        plot(Z(:,2), Z(:,3), 'b.');
        ylabel('z');
        xlabel('y');
    elseif (isImportant(1) == isImportant(3))
        plot(Z(:,1), Z(:,3), 'b.');
        xlabel('x');
        ylabel('z');
    else
        plot(Z(:,1), Z(:,2), 'b.');
        xlabel('x');
        ylabel('y');
    end;
elseif (sum(isImportant) == 1)
    if (isImportant(1) == 1)
        plot(Z(:,1), zeros(1,n), 'm*');
         xlabel('x');
    elseif (isImportant(2) == 1)
        plot(Z(:,2), zeros(1,n), 'm*');
        xlabel('y');
    else
        plot(Z(:,3), zeros(1,n), 'm*');
        xlabel('z');
    end;
end;
grid on;
