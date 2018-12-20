a1c = [9.3;8.2;8;7.3;7.4;6.9;7.3;6.4;6.7;6.7;7.1;6.7];
time_a1c = [
    datetime(2016,03,07,00,00,00);
    datetime(2016,06,20,00,00,00);
    datetime(2016,10,31,00,00,00);
    datetime(2016,12,19,00,00,00);
    datetime(2017,05,04,00,00,00);
    datetime(2017,09,19,00,00,00);
    datetime(2017,10,05,00,00,00);
    datetime(2017,11,22,00,00,00);
    datetime(2017,12,13,00,00,00);
    datetime(2017,12,29,00,00,00);
    datetime(2018,01,30,00,00,00);
    datetime(2018,04,09,00,00,00);
    ];


a1c = a1c - 1.627;
% a1c_train = a1c(7:12);
% time_a1c_train = time_a1c(7:12);
% a1c_test = a1c(1:6);
% time_a1c_test = time_a1c(1:6);
a1c_train = a1c(8);
time_a1c_train = time_a1c(8);

for supp = 108
% Suppose distribution support is 120 days and does not change
% 17 days as an interval; 7 values specified
% m = length(a1c); n = floor(120/m); 
% supp = 101;
m = length(a1c_train); inte = 1; n = floor(supp/inte); 
G_train = zeros(m,n);
for a1c_idx = 1:m
low_bd = time_a1c_train(a1c_idx)-days(supp);
high_bd = time_a1c_train(a1c_idx);
[ tArray, gval_mean ] = mean_glucose( dtime, gval, inte, low_bd, high_bd);
% postprocessing: replace NaN with mean value
for i = 1:length(gval_mean)
    if isnan(gval_mean(i))
        gval_mean(i) = nanmean(gval_mean);
    end
end
% length(gval_mean)
G_train(a1c_idx,:) = flip(gval_mean);
end



v = 20/(250*120)*(inte); % glycation speed every inte days
mu = 1;
G_train = v*G_train;

%(G'*G+mu*eye(size(G,2)))\(G'*a1c)
A = eye(n) + diag(-1*ones(n-1,1),1);

% b = supp:-inte:1; b = flip(b)';
b = 1:inte:supp; b = b';
if length(b) ~= n
    b = b(1:end-1);
end
% k = 1/100/supp;
k = 1e-10;

% for k = 
% cvx_precision low
cvx_begin quiet
    variable f(n+1) nonnegative
%     variable k nonnegative
    minimize norm(G_train*f(1:n) - a1c_train,2) + 1e1*norm(A*f(1:n)-exp(-k*b)*k/(1-exp(-k*supp))*inte,2)
    f(n+1) == 0
    f(1) == 1.0
cvx_end

train_loss(supp) = norm(G_train*f(1:n) - a1c_train,2) + norm(A*f(1:n)-exp(-k*b)*k/(1-exp(-k*supp))*inte,2);
train_a1c_loss(supp) = norm(G_train*f(1:n) - a1c_train,2);

%%%%%%%%%%%%
% subplot(1,3,1)
plot(A*f(1:n))
hold on
plot(exp(-k*b)*k/(1-exp(-k*supp))*inte)
title('C=10')
legend('fitted density','prior density')
return
%%%%%%%%%%%%

%%%%%%%%%%%%%%%% test
%{
m = length(a1c_test); 
G_test = zeros(m,n);
for a1c_idx = 1:m
low_bd = time_a1c_test(a1c_idx)-days(supp);
high_bd = time_a1c_test(a1c_idx);
[ tArray, gval_mean ] = mean_glucose( dtime, gval, inte, low_bd, high_bd);
% postprocessing: replace NaN with mean value
for i = 1:length(gval_mean)
    if isnan(gval_mean(i))
        gval_mean(i) = nanmean(gval_mean);
    end
end
% length(gval_mean)
G_test(a1c_idx,:) = flip(gval_mean);
end

G_test = v*G_test;

test_a1c_loss(supp) = norm(G_test*f(1:n) - a1c_test,2);
%}
end

% train_loss = train_loss(90:130);
% test_a1c_loss = test_a1c_loss(90:130);
train_a1c_loss(supp)
train_loss = train_a1c_loss(90:130);
figure();
plot(90:130,train_loss,'-k')
hold on
plot(90:130,test_a1c_loss,'-.k')
title('training loss and test loss stratified 1');
xlabel('suppport size (days)');
ylabel('loss');
legend('training loss', 'test loss');
