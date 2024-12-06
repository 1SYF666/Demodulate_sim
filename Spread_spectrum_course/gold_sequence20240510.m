clc;close all;clear;
%%
% f(x) = x^6 + x + 1;               1 0 0 0 0 1 1 --->103
% f(x) = x^6 + x^5 + x^2 + x + 1;   1 1 0 0 1 1 1 --->147
%% m1
number_octal1 = 103;
[Binary1,Binary_temp1] = octal_to_binary(number_octal1,8,2);
[m1_out,m1,length_m] = m_generate(Binary_temp1);
%% m2
number_octal2 = 147;

[Binary2,Binary_temp2] = octal_to_binary(number_octal2,8,2);

[m2_out,m2,~] = m_generate(Binary_temp2);

%% 自相关函数
m1 = 1 - msequence*2;

m_temp = m1;

for ii = 1 :length(m1)*10 
    m_temp_register = m_temp(1);
    m_temp(1:end-1) = m_temp(2:end);       % 左移
    m_temp(end) = m_temp_register;
    out(ii) = sum(m1 .* m_temp)/length_m;
end
figure;plot(out);axis([1 length(m1)*3 -0.5 2]);grid on;
title("m1自相关");



