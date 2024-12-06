clc;close all;clear;
%% 参考资料
% C.S0002-0_v1.0中3-17，figure 3.1.3.1.1.1-19
% I and Q Mapping (OTD Mode) for Spreading Rate 1 
% 增加了多径，---20240621

%% 参数设置
Hz = 1 ;            % 单位1Hz
time  = 1;          % 仿真长度
sr = 5000 * Hz; 
ml=2;    
br=sr .* ml; 
nd = sr * time; 
SNR = 5;
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

%% 信道 -- 多径干扰
% 产生多径系数
L = 3;
h_temp = rand(1,L)+1i*rand(1,L);
sum_h = sum(abs(h_temp));
h = h_temp/sum_h;

% 产生时延
t = randi([1 10],1,L-1);

% 产生多径信号
[t_max,index_max] = max(t);
s = scramble_s;
s1_temp = [s zeros(1,t_max)];  
s2_temp = [zeros(1,t(1)) s zeros(1,t_max-t(1))];  
s3_temp = [zeros(1,t(2)) s zeros(1,t_max-t(2))];  

send_signal = h(1).*s1_temp + h(2).*s2_temp + h(3).*s3_temp;

% awgn信道
send_signal = awgn(send_signal,SNR,'measured');

%% rake接收机

send_signal_temp(1,:) = send_signal;
send_signal_temp(2,:) = [send_signal(t(1)+1 : end) zeros(1,t(1))];
send_signal_temp(3,:) = [send_signal(t(2)+1 : end) zeros(1,t(2))];

for k = 1 : L
    % 解扰
    de_scramble = m1 - 1i * m2;
    temp  = [];
    for i = 1 : len
        spread_y_temp = send_signal_temp(k , 1+ (i-1)*N2 : i*N2);
        temp = [temp  spread_y_temp.* de_scramble];
    end
    de_scramble_s(k,:) = temp/2;     % 幅值归一化
    
    % 解扩
    temp_despread = [];
    nd1 =  floor(length(temp)/N1);  % 解扩符号数
    for i = 1 : nd1
        temp_despread = [temp_despread sum(walsh.*de_scramble_s( k , 1 +(i-1)* N1 : i * N1) )/N1];
    end
    de_spread_y(k,:) = temp_despread ;
    clear temp;
    clear temp_despread;
end

% figure;plot(real(de_spread_y));hold on;plot(real(Y)*0.85);

% 多径求和
for k = 1 : L
    r_temp(k,:) =  de_spread_y(k,:).*conj(h(k));
end
r =  sum(r_temp);
%% 解调和判决
r_I = real(r);
r_Q = imag(r);
de_data = [];

for k = 1 : length(r)
    if r_I(k)>0 && r_Q(k)>0
        de_data = [de_data 0 0];
    elseif r_I(k)<0 && r_Q(k)>0
        de_data = [de_data 0 1];
    elseif r_I(k)<0 && r_Q(k)<0
        de_data = [de_data 1 1];
    elseif r_I(k)>0 && r_Q(k)<0
        de_data = [de_data 1 0];
    end
end

%% 误码率计算
errorbit = sum(data(1:length(de_data))~=de_data);
errorrate = errorbit/length(de_data);
fprintf("误码率是：     %f\n",errorrate);
