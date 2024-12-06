%% 频偏估计代码
clc;close all;clear;
%% 参数设置
j = sqrt(-1);
SNR = 10;      % 信噪比
fs = 20e3;
Ts = 1/fs;
N = 31;        % m序列长度
L = 2;         % L段m序列
fc = 100;
%% m序列
number_octal1 = 45;
[Binary1,Binary_temp1] = octal_to_binary(number_octal1,8,2); % 八进制到二进制转换
[m1,~,length_m] = m_generate(Binary_temp1);                  % m序列产生       

%% 导频序列
pilot = [];
for i = 1 : L
   pilot = [pilot m1];
end
%% 信道
r = awgn(pilot,SNR,'measured');
% r = pilot;
r_carrier = r.* exp(j*2*pi*fc*Ts*(1:length(r)));

%% 估计频偏
temp_sum = 0;
for i =  1: N
    temp1 =  r_carrier(i+N*(L-2));
    temp2 =  r_carrier(i+N*(L-1));
%     temp = temp1*conj(temp2);
    temp = temp2*conj(temp1);
    temp_sum = temp_sum + temp;
    theta(i) = angle(temp);
end
fc_est1 =  mean(theta)/(N*Ts*2*pi);    % 方法一
fc_est2 = angle(temp_sum)/(N*Ts*2*pi); % 方法二
fprintf("方法一估计频偏：%f\n",fc_est1);
fprintf("方法二估计频偏：%f\n",fc_est2);