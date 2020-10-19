function [] = foo1(p, mu1, mu2, Sigma, n1, n2, n)


%Итак, приступим к построению обучающих выборок

format compact


q1 = n1 / (n1 + n2);
q2 = n2 / (n1 + n2);

%генерируем TS1 и TS2 (нормальное трехмерное распределение)
TS1 = mvnrnd(mu1, Sigma, n1);
TS2 = mvnrnd(mu2, Sigma, n2);

% Построение обучающих выорок завершено
%Создадим тестовую выборку размера 2n
%n = 512;
N = 2 * n;

selection1 = mvnrnd(mu1, Sigma, n);
selection2 = mvnrnd(mu2, Sigma, n);

test_selection = cat(1, selection1, selection2);

%Найдем оценки параметров распределений
m_estimate_1 = [0 0 0];
m_estimate_2 = [0 0 0];
for j = 1:p
    m_estimate_1(j) = mean(TS1(:, j));
    m_estimate_2(j) = mean(TS2(:, j));
end

% найдём матрицу S, которая является оценкой матрицы ковариаций Sigma
S1 = zeros(p);
S2 = zeros(p);

for i = 1:p
    for j = 1:p
        S1(i, j) = sum((TS1(:, i) - m_estimate_1(i)) .* (TS1(:, j) - m_estimate_1(j))) / (n1 - 1);
    end
end
%-------------------------------------------------------------------------------------------------------------------
for i = 1:p
    for j = 1:p
        S2(i, j) = sum((TS2(:, i) - m_estimate_2(i)) .* (TS2(:, j) - m_estimate_2(j))) / (n2 - 1);
    end
end

S = ((n1 - 1) * S1 + (n2 - 1) * S2) / (n1 + n2 - 2);

% Найдём оценку вектора альфа
a = S \ ((m_estimate_1 - m_estimate_2)');

z1 = TS1 * a;
z2 = TS2 * a;

z_estimate_1 = mean(z1);
z_estimate_2 = mean(z2);

s2_z = 0;
for i = 1:p
    for j = 1:p
        s2_z = s2_z + a(i) * S(i, j) * a(j);
    end
end

%классифицируем тестовую выборку
result = 1:N;
rigth_side = (z_estimate_1 + z_estimate_2) / 2 + log(q2 / q1);
for i = 1:N
    if (test_selection(i, :) * a >= rigth_side)
        result(i) = 1;
    else
        result(i) = 2;
    end
end

disp('Biased estimate of Mahalanobis distance:');
biased_estimate_of_Mahalanobis_distance = (z_estimate_1 - z_estimate_2).^2 / s2_z;
disp(biased_estimate_of_Mahalanobis_distance);
disp('Unbiased estimate of the Mahalanobis distance:');
unbiased_estimate_of_Mahalanobis_distance = (n1 + n2 - p - 3) * biased_estimate_of_Mahalanobis_distance / (n1 + n2 - 2) - p * ((1 / n1) + (1 / n2));
disp(unbiased_estimate_of_Mahalanobis_distance);

%результаты классификации тестовой выборки
four_field_contingency_table = zeros(2);
true1 = 0;
true2 = 0;
false1 = 0;
false2 = 0;
for i = 1:n
    if (result(i) == 1)
        true1 = true1 + 1;
    else
        false2 = false2 + 1;
    end
end
for i = n+1:N
     if (result(i) == 1)
        false1 = false1 + 1;
    else % result == 2
        true2 = true2 + 1;
    end 
end
four_field_contingency_table(1, 1) = true1;
four_field_contingency_table(1, 2) = false2;
four_field_contingency_table(2, 1) = false1;
four_field_contingency_table(2, 2) = true2;

disp('Classification results of test_set:');
disp(four_field_contingency_table)

disp('Probability error:')
disp('P(2|1):')
P_roof21 = false2 / (true1 + false2);
disp(P_roof21);
disp('P(1|2):')
P_roof12 = false1 / (true2 + false1);
disp(P_roof12);

%классифицируем TS1 и TS2
OB = cat(1, TS1, TS2);
result = 1:(n1 + n2);
for i = 1 : (n1 + n2)
    if ( OB(i, :) * a >= rigth_side)
        result(i) = 1;
    else
        result(i) = 2;
    end
end

%результаты классификации TS1 и TS2
four_field_contingency_table = zeros(2, 2);
true1 = 0;
true2 = 0;
false1 = 0;
false2 = 0;
for i = 1 : n1
    if (result(i) == 1)
        true1 = true1 + 1;
    else % result(i) == 2
        false2 = false2 + 1;
    end
end
for i = (n1 + 1) : (n1 + n2)
    if (result(i) == 1)
        false1 = false1 + 1;
    else % result(i) == 2
        true2 = true2 + 1;
    end 
end
%true1 + true2 + false1 + false2 == (n1 + n2);
four_field_contingency_table(1, 1) = true1;
four_field_contingency_table(1, 2) = false2;
four_field_contingency_table(2, 1) = false1;
four_field_contingency_table(2, 2) = true2;

disp('Classification result of TS1 and TS2:');
disp(four_field_contingency_table);

%оценим вероятность ошибочной классификации
K = log(q2 / q1);
P_roof21 = normcdf( (( K - unbiased_estimate_of_Mahalanobis_distance / 2 ) / sqrt(unbiased_estimate_of_Mahalanobis_distance)) , 0 , 1 );
P_roof12 = normcdf( ((-K - unbiased_estimate_of_Mahalanobis_distance / 2) / sqrt(unbiased_estimate_of_Mahalanobis_distance)), 0 , 1);
disp('Probability of error:')
disp('P(2|1):')
disp(P_roof21);
disp('P(1|2):')
disp(P_roof12);
