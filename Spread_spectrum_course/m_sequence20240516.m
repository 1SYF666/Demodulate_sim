clc;close all;clear;
%% m1
number_octal1 = 45;
[Binary1,Binary_temp1] = octal_to_binary(number_octal1,8,2);
[m1,~,length_m] = m_generate(Binary_temp1);
%% m2
number_octal2 = 57;
[Binary2,Binary_temp2] = octal_to_binary(number_octal2,8,2);
[m2,~,~] = m_generate(Binary_temp2);
%% 自相关
m_temp = m1;
for ii = 1 :length(m_temp)*10 
    m_temp_register = m_temp(1);
    m_temp(1:end-1) = m_temp(2:end);       % 左移
    m_temp(end) = m_temp_register;
    out(ii) = sum(m1 .* m_temp)/length_m;
end

figure;plot(out);axis([1 length(m1)*3 -0.2 1.2]);grid on;

%% 互相关
m_temp = m2;
for ii = 1 :length(m_temp)*10
    m_temp_register = m_temp(1);
    m_temp(1:end-1) = m_temp(2:end);       % 左移
    m_temp(end) = m_temp_register;
    out2(ii) = sum(m1 .* m_temp)/length_m;
end
figure;plot(out2);
