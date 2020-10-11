K = 6;
%disp('foo1 is starting');
%foo1(3, [0, 0, 0], [0.3, 0.8, 1], eye(3), );
%disp('foo1 is finishing');
%disp('------------------------------------------------------------------------')
%disp('------------------------------------------------------------------------')
%disp('foo2 is starting');
I = zeros(1, K);
D = zeros(1, K);
D_H = zeros(1, K);
FCT_Set = zeros(2);
FCT_Set(:, :, K) = zeros(2);
FCT_TS = zeros(2);
FCT_TS(:, :, K) = zeros(2);
P_21 = zeros(3, K);
P_12 = zeros(3, K);
for i = 1:K
%i = 1;
%    disp('50i == ');
%    disp(50*i);
    [D(i), D_H(i), FCT_Set(:,:,i), FCT_TS(:,:,i), P_21(:, i), P_12(:,i)] = foo2(50*i);
    I(i) = 50*i;
end;
TableForPrint = cat(2, I', D', D_H', (D - D_H)');
TableForPrint_2 = cat(2, I', P_21', P_12');
disp('-------------------------------------')
disp(TableForPrint);
disp('---')
disp(TableForPrint_2);
disp('-------------------------------------')
%for i = 1:K
%    disp('FCT_Set');
%    disp(FCT_Set(:,:,i));
%    disp('FCT_TS');
%    disp(FCT_TS(:,:,i));
%end;
disp('foo2 is finishing');