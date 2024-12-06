clc;
close all;
clear;
%%
mu = 3.8;
MaxCycle = 50;
x = zeros(1,MaxCycle);
x(1) = 0.5;  %  初始值
i = 1;
while i <= MaxCycle - 1
    x(i + 1) = x(i) * mu * (1 - x(i));
    i = i + 1;
end
% 绘图
for ii = 1 :length(x) 
    xcorrout2(ii) = sum(x.*circshift(x,ii-1))/length(x);
end
figure;plot(xcorrout2);grid on;
title("自相关");

% 阈值
x2 = (xcorrout2>0.4);
for ii = 1 :length(x2) 
    xcorrout3(ii) = sum(x2.*circshift(x2,ii-1))/length(x2);
end
figure;plot(xcorrout3);grid on;
title("自相关");


