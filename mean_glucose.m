function [ tArray, gval_mean ] = mean_glucose( dtime, gval, int, low_bd, high_bd )
% mean glucose every interval
if nargin <= 3
low_bd = min(dtime);
high_bd = max(dtime);
end
tArray = low_bd:days(int):high_bd;
% size(tArray)
gval_mean = zeros(size(tArray,2)-1,1);
for i=1:length(tArray)-1
    t_low = datetime(tArray(i));
    t_high = datetime(tArray(i+1));
    idx = find((dtime>=t_low)&(dtime<t_high));
%     dtime_int = dtime(idx);
    gval_int = gval(idx);
    gval_mean(i) = nanmean(gval_int);
%     if isnan(gval_mean(i))
%         gval_mean(i) = gval_mean(i-1);
%     end
end
tArray = tArray(1:end-1);
end

