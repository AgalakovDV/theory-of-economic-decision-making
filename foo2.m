%number - размер обучающей выборки. 
function [D, D_H, FCT_Set, FCT_TS, P_21, P_12] = foo2(number)

format compact
FCT_Set = zeros(2);
FCT_TS = zeros(2);
P_21 = zeros(1,3);
P_12 = zeros(1,3);
filename = 'D:/7 семестр/видеолекции/Теория Принятия Экономических Решений/лабораторные работы/лаба1/german/german.data-numeric';
delimiterIn = ' ';
headerlinesIn = 0;
A = importdata(filename, delimiterIn, headerlinesIn);
p = 24;

%disp('A.size')
%disp(size(A))

%построение ОВ
TS11 = A(A(:, 25) == 1, :);
TS22 = A(A(:, 25) == 2, :);

TS1 = TS11(1:number, :);
TS2 = TS22(1:number, :);

%disp('size helper');
%disp(size(TS11));
%disp('number == ');
%disp(number);
%disp('size');
%disp(size(TS1));

test_selection = cat(1, TS11(number+1:size(TS11, 1), :), TS22(number + 1:size(TS22, 1), :));

n1 = size(TS1, 1);
n2 = size(TS2, 1);
N = size(test_selection, 1);

%q1 = n1 / (n1 + n2);
%q2 = n2 / (n1 + n2);
K = log(n2 / n1);%log(q2/q1)

m_estimate_1 = zeros(1,p);
m_estimate_2 = zeros(1,p);
for j = 1:p
    m_estimate_1(j) = sum(TS1(:, j)) / n1;
    m_estimate_2(j) = sum(TS2(:, j)) / n2;
end

S1 = zeros(p);
S2 = zeros(p);

for l = 1:p
    for j = 1:p
        S1(l, j) = sum((TS1(:, l) - m_estimate_1(l)) .* (TS1(:, j) - m_estimate_1(j))) / (n1 - 1);
    end
end

for l = 1:p
    for j = 1:p
        S2(l, j) = sum((TS2(:, l) - m_estimate_2(l)) .* (TS2(:, j) - m_estimate_2(j))) / (n2 - 1);
    end
end

S = ((n1 - 1) * S1 + (n2 - 1) * S2) / (n1 + n2 - 2);

%заменяем альфа на а
a = S \ ((m_estimate_1 - m_estimate_2)');

%zs1 = sum(TS1) / n1;
%zs2 = sum(TS2) / n2;
%disp('zs1');
%disp(zs1');
%disp('zs2'); 
%disp(zs2');


z1 = TS1(:, 1:p) * a;
z2 = TS2(:, 1:p) * a;

%disp(z1);

z_roof1 = sum(z1) / n1;
z_roof2 = sum(z2) / n2;

%disp('Z_roof');
%disp(z_roof1);
%disp(z_roof2);
%disp(n1);
%disp(n2);

s2_z = 0;
for l = 1:p
    for j = 1:p
        s2_z = s2_z + a(l) * S(l, j) * a(j);
    end
end

%классифицируем элементы тестовой выборки
class = zeros(1,N);
C = (z_roof1 + z_roof2) / 2 + K;
for i = 1:N
    if (test_selection(i, 1:p) * a >= C)
        class(i) = 1;
    else
        class(i) = 2;
    end
end

D = (z_roof1 - z_roof2).^2 / s2_z;
D_H = (n1 + n2 - p - 3) * D(1) /(n1 + n2 - 2) - p * ((1 / n1) + (1 / n2));

%результаты классификации тестовой выборки
for i = 1:N
    if (class(i) == 1 && test_selection(i, 25) == 1)
        FCT_Set(1, 1) = FCT_Set(1, 1) + 1;
    end
     if (class(i) == 2 && test_selection(i, 25) == 1)
        FCT_Set(1, 2) = FCT_Set(1, 2) + 1;
    end
    if (class(i) == 1 && test_selection(i, 25) == 2)
        FCT_Set(2, 1) = FCT_Set(2, 1) + 1;
    end
    if (class(i) == 2 && test_selection(i, 25) == 2)
        FCT_Set(2, 2) = FCT_Set(2, 2) + 1;
    end 
end

%классифицируем ОВ1 и ОВ2
TS = cat(1, TS1, TS2);
class = zeros(1,n1+n2);
for i = 1 : (n1 + n2)
    if ( TS(i, 1:p) * a >= C)
        class(i) = 1;
    else
        class(i) = 2;
    end
end

%результаты классификации ОВ1 и ОВ2
for i = 1 : (n1 + n2)
    if (class(i) == 1 && TS(i, 25) == 1)
        FCT_TS(1, 1) = FCT_TS(1, 1) + 1;
    end
    if (class(i) == 2 && TS(i, 25) == 1)
        FCT_TS(1, 2) = FCT_TS(1, 2) + 1;
    end
    if (class(i) == 1 && TS(i, 25) == 2)
        FCT_TS(2, 1) = FCT_TS(2, 1) + 1;
    end
    if (class(i) == 2 && TS(i, 25) == 2)
        FCT_TS(2, 2) = FCT_TS(2, 2) + 1;
    end 
end


%оценим вероятность ошибочной классификации

P_21(1) = FCT_Set(2, 1) / (FCT_Set(1, 1) + FCT_Set(2, 1));
P_12(1) = FCT_Set(1, 2) / (FCT_Set(2, 2) + FCT_Set(1, 2));
P_21(2) = normcdf( (( K - D_H / 2 ) / sqrt(D_H)) , 0 , 1 );
P_12(2) = normcdf( ((-K - D_H / 2) / sqrt(D_H)), 0 , 1);
P_21(3) = FCT_TS(2, 1) / (FCT_TS(1, 1) + FCT_TS(2, 1));
P_12(3) = FCT_TS(1, 2) / (FCT_TS(2, 2) + FCT_TS(1, 2));


%disp('Biased estimate of Mahalanobis distance:');
%disp(D);
%disp('Unbiased estimate of the Mahalanobis distance:');
%disp(D_H);
%disp('Classification classs of test_set:');
%disp(FCT_Set);
%disp('Probability error:')
%disp('P(2|1):')
%disp(P_21(1));
%disp('P(1|2):')
%disp(P_12(1));
%disp('Classification class of TS1 and TS2:');
%disp(FCT_TS);
%disp('Difference probability of error:');
%disp('P(2|1):');
%disp(P_21(2) - P_21(3));
%disp('P(1|2):');
%disp(P_12(2));
%disp(P_12(3));