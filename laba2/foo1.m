function [] = foo1(lvl)
 
% f = fopen(filename_, 'rt');
% while (feof ~= true)
%     datas = textscan(f, '%f %s %f %s %f %s');
% end;
% fclose(f);
% celldisp(datas);

% fileID = fopen(filename,'r');
% dataArray = fscanf(inFile,'%f %f %f %f %s',[2 Inf])';
% fclose(fileID);


format compact;
filename = 'D:/7 семестр/видеолекции/Теория Принятия Экономических Решений/лабораторные работы/лаба2/iris/iris2.data';
f1 = 'D:/7 семестр/видеолекции/Теория Принятия Экономических Решений/лабораторные работы/лаба2/iris/ir1.txt';
f2 = 'D:/7 семестр/видеолекции/Теория Принятия Экономических Решений/лабораторные работы/лаба2/iris/ir2.txt';
f3 = 'D:/7 семестр/видеолекции/Теория Принятия Экономических Решений/лабораторные работы/лаба2/iris/ir3.txt';
delimiterIn = ',';
headerlinesIn = 0;
set = importdata(filename, delimiterIn, headerlinesIn);
p = 4;
set1 = importdata(f1, delimiterIn, headerlinesIn); 
set2 = importdata(f2, delimiterIn, headerlinesIn);
set3 = importdata(f3, delimiterIn, headerlinesIn);

%центрируем матрицу set
mean_vec = mean(set);
for i = 1:p
    set(:,i) = set(:,i) - mean_vec(i);
end;

Sigma = cov(set);

vec = (size(set));
n = vec(1);
vec = (size(set1));
n1 = vec(1);
vec = (size(set2));
n2 = vec(1);
vec = (size(set3));
n3 = vec(1);
% Приступим к переходу от x к z
Z = zeros(n, p);
Z1 = zeros(n1, p);
Z2 = zeros(n2, p);
Z3 = zeros(n3, p);
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
        Z(k, i) = set(k,:) * L(:, ind(p + 1 - i));%
    end;
end;

for k = 1:n1
    for i = 1:p
        Z1(k, i) = set1(k,:) * L(:, ind(p + 1 - i));%
    end;
end;

for k = 1:n2
    for i = 1:p
        Z2(k, i) = set2(k,:) * L(:, ind(p + 1 - i));%
    end;
end;

for k = 1:n3
    for i = 1:p
        Z3(k, i) = set3(k,:) * L(:, ind(p + 1 - i));%
    end;
end;

%disp('Set = ');
%disp(set);
%disp('Z = ');
%disp(Z);
%disp(L);
 % Найдём матрицу ковариаций по z
 Sigma_z2 = eye(p);
 for i = 1:p
     Sigma_z2(i,i) = r_sort(p + 1 - i);
 end;

 
Sigma_proc = zeros(1,p);
sum0 = 0;
for i = 1:p
    Sigma_proc(i) = sum0 + Sigma_z2(i,i) / sum(r_sort);
    sum0 = Sigma_proc(i);
end;

disp('Sigma_proc');
disp(Sigma_proc');

% Порог lvl%
isImportant = zeros(1,p);
for i = 1:p
    isImportant(i) = 1;
    if (Sigma_proc(i) > lvl)
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

figure('Name','Measured Data');
hold on;
title('Z_3');
if (sum(isImportant) == 2)
    plot(Z1(:,1), Z1(:,2), 'b.');
    plot(Z2(:,1), Z2(:,2), 'm.');
    plot(Z3(:,1), Z3(:,2), 'g.');
    xlabel('x');
    ylabel('y');
elseif (sum(isImportant) == 1)
    plot(Z(:,1), zeros(1,n), 'm*');
    xlabel('x');
elseif (sum(isImportant) == 3)
    plot3(Z1(:,1), Z1(:,2), Z1(:,3), 'b.');
    plot3(Z2(:,1), Z2(:,2), Z2(:,3), 'm.');
    plot3(Z3(:,1), Z3(:,2), Z3(:,3), 'g.');
    view(50, 40);
    xlabel('x');
    ylabel('y');
    zlabel('z');
end;
grid on;
