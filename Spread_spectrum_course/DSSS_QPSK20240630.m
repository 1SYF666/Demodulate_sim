clc;close all;clear;
%% 参考资料
% C.S0002-0_v1.0中3-17，figure 3.1.3.1.1.1-19
% I and Q Mapping (OTD Mode) for Spreading Rate 1 
%% notice
% 增加了多径，---20240621
% 修改了解扰解扩数据长度，按照多少帧截取处理  ---20240630
% 增加了了低通滤波  ---20240630

%% 参数设置
Hz = 1 ;            % 单位1Hz
time  = 1;          % 仿真长度
sr = 5000 * Hz; 
ml=2;    
br=sr .* ml; 
nd = sr * time;
% SNR = 5;
EbN0 = 0 : 7;
% 20240630
walsh_len = 64;             % 扩频码长度
r = 10;
scramble_len = 2^r - 1;     % 扰码长度 
Rc = sr * walsh_len;        % 码片速率
sps = 5;                    % 每个码片采样四个点
fs = sps * Rc;              % 采样率

%% 基带信息
data=rand(1,nd*ml)>0.5;
% QPSK 调制 
Y = zeros(1,nd);
for i = 1 : nd
    if data(1+ ml*(i-1) : ml*i) == [0 0]
        Y(i) = 1 + 1i * 1 ;
    elseif data(1+ ml*(i-1) : ml*i) == [0 1]
        Y(i) = -1 + 1i * 1 ;
    elseif data(1+ ml*(i-1) : ml*i) == [1 1]
        Y(i) = -1 - 1i * 1 ;
    elseif data(1+ ml*(i-1) : ml*i) == [1 0]
        Y(i) = 1 - 1i * 1 ;
    end
end

Y = Y / sqrt(2);  % 功率归一化

% figure;scatter(real(Y),imag(Y));title("QPSK调制后星座图");

%% 扩频
% wash -- 64阶
hadamard = h_generate(walsh_len);
walsh = hadamard(15,:);
spread_y = [];
for i = 1 : nd
    spread_y = [spread_y walsh*Y(i)];
end

%% 加扰
number_octal1 = 2011; % 10阶
[Binary1,Binary_temp1] = octal_to_binary(number_octal1,8,2);
[m1,~,length_m] = m_generate(Binary_temp1);  % m1

number_octal2 = 1881; % 10阶
[Binary2,Binary_temp2] = octal_to_binary(number_octal2,8,2);
[m2,~,length_m2] = m_generate(Binary_temp2); % m1
scramble = m1 + 1i * m2;

len = floor(length(spread_y)/scramble_len); % 加扰信号片数
scramble_y = [];
for i = 1 : len
    spread_y_temp = spread_y(1+ (i-1)*scramble_len : i*scramble_len);
    scramble_y = [scramble_y  spread_y_temp.* scramble];
end

% figure;scatter(real(scramble_s),imag(scramble_s));title("加扩加扰后星座图");

%% 成型滤波
s_real = real(scramble_y);
s_imge = imag(scramble_y);

rolloff_factor = 0.5;       % 滚降因子
rcos_fir = rcosdesign(rolloff_factor,2*sps,sps); % 默认是根升余弦滤波器,'sqrt'

%%% 插值
% up_sps_ds = zeros(1,length(ds)*sps);
for k=1:length(s_real)
    up_sps_ds_real(1+sps*(k-1)) = s_real(k);
    up_sps_ds_real(2+sps*(k-1):sps*k) = zeros(1,sps-1);
end

for k=1:length(s_imge)
    up_sps_ds_imag(1+sps*(k-1)) = s_imge(k);
    up_sps_ds_imag(2+sps*(k-1):sps*k) = zeros(1,sps-1);
end

%%% 滚降滤波
rcos_ds_real = filter(rcos_fir,1,up_sps_ds_real);
rcos_ds_imag = filter(rcos_fir,1,up_sps_ds_imag);

scramble_y = rcos_ds_real + sqrt(-1)*rcos_ds_imag;

%% 信道 -- 多径干扰
% 产生多径系数
L = 3;
h_temp = rand(1,L)+1i*rand(1,L);
sum_h = sum(abs(h_temp));
h = h_temp/sum_h;  % 多径系数归一化操作

% 产生时延
t = randi([1 10],1,L-1);

% 产生多径信号
[t_max,index_max] = max(t);
s = scramble_y;
s1_temp = [s zeros(1,t_max)];  
s2_temp = [zeros(1,t(1)) s zeros(1,t_max-t(1))];  
s3_temp = [zeros(1,t(2)) s zeros(1,t_max-t(2))];  

send_signal = h(1).*s1_temp + h(2).*s2_temp + h(3).*s3_temp;

% awgn信道
Ebn0(1) = EbN0(1)-log(sps)/log(10);
send_signal = awgn(send_signal,Ebn0(1),'measured');

%% 匹配滤波



%% rake接收机

frame_symbol = 3200 ;                                 % 假设一帧长度为framelen个符号
spread_frame_len = frame_symbol*walsh_len;
send_signal_frame = send_signal(1 : spread_frame_len); % 截取一帧长度数据进行解扰解扩
len = floor(length(send_signal_frame)/scramble_len);   % 一帧里面有多少片扰码

for k = 1 : L

    de_scramble = m1 - 1i * m2;                 % 解扰
    temp  = [];
    
    if k == 1                                   % 去除多径延迟
        send_signal_temp = send_signal_frame;
    elseif k==2
        send_signal_temp = send_signal_frame(t(1)+1 : end);
    else
        send_signal_temp = send_signal_frame(t(2)+1 : end);
    end
    
    for i = 1 : len
        spread_y_temp = send_signal_temp(1+ (i-1)*scramble_len : i*scramble_len);
        temp = [temp  spread_y_temp.* de_scramble];
    end
    de_scramble_s(k,:) = temp/2;                % 幅值归一化
    
    
    temp_despread = [];                         % 解扩
    nd1 =  floor(length(temp)/walsh_len);       % 解扩符号数
    for i = 1 : nd1
        temp_despread = [temp_despread sum(walsh.*de_scramble_s( k , 1 +(i-1)* walsh_len : i * walsh_len) )/walsh_len];
    end
    de_spread_y(k,:) = temp_despread ;
    clear temp;
    clear temp_despread;
end % k = 1 : L

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
