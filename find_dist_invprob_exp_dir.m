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
a1c_train = a1c(1:6);
time_a1c_train = time_a1c(1:6);
% a1c_test = a1c(7:12);
% time_a1c_test = time_a1c(7:12);
a1c_train = a1c(10:11);
time_a1c_train = time_a1c(10:11);

TOTAL_LOSS = [];
Ks = [];

for supp = 100:130
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
b = 0:inte:supp; b = b';
if length(b) ~= n
    b = b(1:end-1);
end
% k = 1/100/supp;
k = 1e-2;

dfdk_com = @(t,k) (supp*exp(-k*supp)*(exp(-k*supp) - exp(-k*t)))./(exp(-k*supp) - 1)^2 - (supp*exp(-k*supp) - t.*exp(-k*t))./(exp(-k*supp) - 1);
% dfdk = deriv_k(dfdk_com, b, k);

f = @(k) (exp(-k*b) - exp(-k*supp)) / (1-exp(-k*supp));
% grad = (G_train*f(k)-a1c_train)'*G_train*dfdk;
df = @(k) k*exp(-k*b)/(1-exp(-k*supp));


loss = [];
i = 1;
while(true)
loss(i) = norm(G_train*f(k) - a1c_train,2);
if loss(i) < 1e-10
    break;
end
grad = grad_k(k,dfdk_com,f,G_train,a1c_train,b);
k = k - (1e-7)*grad;
i = i + 1;
if i > 20000
    break;
end
end

k
loss(end)

TOTAL_LOSS(supp) = loss(end);
Ks(supp) = k;



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

test_loss(supp) = norm(G_test*f(1:n) - a1c_test,2);
%}
end

train_loss = train_loss(90:130);
% test_loss = test_loss(90:130);
figure();
plot(90:130,train_loss,'-')
% hold on
% plot(90:130,test_loss,'-.')
title('training loss on whole data');
xlabel('suppport size (days)');
ylabel('loss');
% legend('training loss', 'test loss');


