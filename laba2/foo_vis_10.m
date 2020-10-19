function [] = foo_vis_10(n)

format compact;
format long;
p = 10;
mu = zeros(1,p);
%          x1  x4 3x1-1 2*(6x1-1)  x2  x2+x1  x3  x2+x3  2x2
m10 =     [ 1   3   2   10  -10     1    2    -3    -2    2;
            2   4   5   22  -22     0    2    -1    -1    0;
            0   5  -1  -2    2      0    0    -9    -9    0;
            3   6   8   34  -34     0    3     9     9    0;
            2   0   5   22  -22    17   19     9    26   34;
           -2  -4  -7  -26   26     6    4     8    14   12;
           -1  12  -4  -14   14     6    5     7    13   12;
           -3  61  -10 -38   38     1   -2     6     7    2;
            0   0  -1  -2    2      1    1     5     6    2;
            1   2   2  10   -10     0    1    17    17    0;
           ];

% m10 = magic(p);
% m10(1,1) = 100;
% m10(10,3) = 17;
% m10(10,10) = 100;
Sigma = cov(m10/10);

% disp('Sigma');
% disp(Sigma);

for i = 1:5
    Sigma(i+5,i+5) = Sigma(i+5,i+5) * 8;
end;

disp('Sigma');
disp(Sigma);

% fname = 'D:/7 семестр/видеолекции/Теория Принятия Экономических Решений/лабораторные работы/лаба2/sigma.xlsx';
% f = fopen(fname, 'wt');
% for i = 1:10
%     for j = 1:10
%         fprintf(f, '%f %s', Sigma(i,j), ' ');
%     end;
%     fprintf(f, '%s', ' \n ');
% end;
% fclose(f);


set = mvnrnd(mu, Sigma, n);

% Приступим к переходу от x к z
Z = zeros(n, p);
% z_i = l_i*x, l_i - с.в. матрицы ковариаций вектора x

% Найдём собственные значения матрицы ковариаций
[L, R] = eig(Sigma);

% Отсортируем с.ч.
% r_sort - массив, отсортированных с.ч.
% ind - массив номеров, в каком порядке они находились.
[r_sort, ind] = sort(diag(R));
for k = 1:n
    for i = 1:p
        Z(k, i) = set(k,:) * L(:, ind(p + 1 - i));
    end;
end;

 % Найдём матрицу ковариаций по z
 Sigma_z2 = eye(p);
 for i = 1:p
     Sigma_z2(i,i) = r_sort(p + 1 - i);
 end;
 
disp('r_sort');
disp(r_sort);
%disp('Sigma_z2')
%disp(Sigma_z2);



% Процентная доля дисперсии на каждую из компонент
% Sigma_proc - массив долей дисперсий для каждого измерения из p штук

Sigma_proc = zeros(1,p);
sum0 = 0;
for i = 1:p
    Sigma_proc(i) = sum0 + Sigma_z2(i,i) / sum(r_sort);
    sum0 = Sigma_proc(i);
end;

disp('Sigma_proc');
disp(Sigma_proc');

% Порог 70%
isImportant = zeros(1,p);
for i = 1:p
    isImportant(i) = 1;
    if (Sigma_proc(i) > 0.7)
        break
    end;
end;


figure('Name','Measured Data');
hold on;
title('Z_2');

if (sum(isImportant) == 2)
    plot(Z(:,1), Z(:,2), 'b.');
    xlabel('x');
    ylabel('y');
elseif (sum(isImportant) == 1)
    plot(Z(:,1), zeros(1,n), 'm*');
    xlabel('x');
elseif (sum(isImportant) == 3)
    plot3(Z(:,1), Z(:,2), Z(:,3), 'gv');
    view(50, 40);
    xlabel('x');
    ylabel('y');
    zlabel('z');
end;
grid on;

