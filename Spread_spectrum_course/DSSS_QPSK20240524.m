clc;close all;clear all;
%% 参考资料
% C.S0002-0_v1.0中3-17，figure 3.1.3.1.1.1-19
% I and Q Mapping (OTD Mode) for Spreading Rate 1 

%% 参数设置
Hz = 1 ;            % 单位1Hz
time  = 1;          % 仿真长度
sr = 50 * Hz; 
ml=2;    
br=sr .* ml; 
nd = sr * time; 

%% 基带信息
data=rand(1,nd*ml)>0.5;
% QPSK 调制 
Y = zeros(1,nd);
for i = 1 : nd
    if data(1+ ml*(i-1) : ml*i) == [0 0]
        Y(i) = 1/sqrt(2) + 1i * 1/sqrt(2);
    elseif data(1+ ml*(i-1) : ml*i) == [0 1]
        Y(i) = -1/sqrt(2) + 1i * 1/sqrt(2);
    elseif data(1+ ml*(i-1) : ml*i) == [1 1]
        Y(i) = -1/sqrt(2) - 1i * 1/sqrt(2);
    elseif data(1+ ml*(i-1) : ml*i) == [1 0]
        Y(i) = 1/sqrt(2) - 1i * 1/sqrt(2);
    end
end
% 星座图
figure;scatter(real(Y),imag(Y));

%% 扩频
% wash -- 64阶
N1 = 64;
hadamard = h_generate(N1);
walsh = hadamard(15,:);
spread_y = [];
for i = 1 : nd
    spread_y = [spread_y walsh*Y(i)];
end

%% 加扰
r = 10;
N2 = 2^r - 1;
number_octal1 = 2011; % 10阶
[Binary1,Binary_temp1] = octal_to_binary(number_octal1,8,2);
[m1,~,length_m] = m_generate(Binary_temp1);  % m1

number_octal2 = 1881; % 10阶
[Binary2,Binary_temp2] = octal_to_binary(number_octal2,8,2);
[m2,~,length_m2] = m_generate(Binary_temp2); % m1
scramble = m1 + 1i * m2;

len = floor(length(spread_y)/N2); % 加扰信号片数
scramble_s = [];
for i = 1 : len
    spread_y_temp = spread_y(1+ (i-1)*N2 : i*N2);
    scramble_s = [scramble_s  spread_y_temp.* scramble];
end

figure;scatter(real(scramble_s),imag(scramble_s));

%% 信道



%% 解扰
de_scramble = m1 - 1i * m2;
de_scramble_s = [];
for i = 1 : len
    spread_y_temp = scramble_s(1+ (i-1)*N2 : i*N2);
    de_scramble_s = [de_scramble_s  spread_y_temp.* de_scramble];
end
de_scramble_s = de_scramble_s/2;     % 幅值归一化

%% 解扩
de_spread_y = [];
nd1 =  floor(length(de_scramble_s)/N1);  % 解扩符号数
for i = 1 : nd1
    de_spread_y = [de_spread_y sum(walsh.*de_scramble_s( 1 +(i-1)* N1 : i * N1) )/N1];
end
figure;plot(real(de_spread_y));hold on;plot(real(Y)*0.85);

%% 解调






