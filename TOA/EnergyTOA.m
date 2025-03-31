clc;clear;
close all;
j = sqrt(-1);
%% 参数设置
K = 1e3;             % 单位 KHz
M = 1e6;             % 单位 MHz
Rb = 1*K;            % 码速率
fs = 4 * Rb;      % 采样率
time = 5;            % 仿真时间
symbolnum = Rb*time; % 码元个数
fc = 1500*M;       
SNR_dB = 15;        % 信噪比(单位 dB);
start=30000;
finish=30000;
T = 1/fs;           % 采样间隔 (s)
starttime=start*T;
t = 0:T:1-T;        % 时间向量 (1s)
% T_b = 0.1;          % 能量累积周期 (s)
% W = round(T_b / T); % 每个能量块的采样点数
W = 32;      %能量块采样点数   
W_slip = 32;   %移动窗长
threadhold1 = 28;   %移动窗检测门槛
toa_est = [];  %起始点
toe_est = [];  %终止点
startpoint = [];
finishpoint = [];
Esum = [];
Ysum = [];
start_counter = 0;  % 起始点计数器
finish_counter = 0; % 终止点计数器

%% BPSK信号生成
% PCM编码
pcm_bits = randi([0 1],1,symbolnum);  %生成长度为symbolnum的随机（randi）二进制序列(pcm_bits)，值为0或1
pcm_symbols = 2 * pcm_bits - 1;   % 将0 --> -1; 1 --> 1   %pcm_symbols,映射后的BPSK序列
% BPSK 调制
sps = fs / Rb;                    % sps每符号采样数
rcos_pcm_s1 = repmat(pcm_symbols,sps,1);     %repmat将每个符号重复sps次
rcos_pcm_s1(2:sps,:) = 0;      %将第二行到第sps行的值设为0
rcos_pcm_s = rcos_pcm_s1(:)';   %将矩阵展平为行向量
rolloff_factor = 0.7;       % 滚降因子，取值范围[0,1]，控制滤波器的带宽和旁瓣衰减
rcos_fir = rcosdesign(rolloff_factor,2*sps,sps); % 默认是根升余弦滤波器,'sqrt'
%2*sps，滤波器的阶数，长度为2*sps+1；rcos_fir，根升余弦滤波器系数
rcos_ds_real = filter(rcos_fir,1,rcos_pcm_s);  %1，滤波器的分母系数；rcos_pcm_real，脉冲成形后的信号

%% 单源信号模型
S1 = rcos_ds_real+j*rcos_ds_real;  %S1为复信号，实部+虚部
%fa = exp(-j*2*pi*fc*d*sin(theta)/c);
%disp(['S1 的长度：', num2str(length(S1))]);  % 打印 S1 的长度
Y = [zeros(1,start) S1  zeros(1,finish) S1 zeros(1,finish) S1 zeros(1,finish)];   %生成发送信号，前后加3W个零（零填充）
%figure;plot(real(Y));   %绘制零填充后Y的实部

%% 信道
SNR1_dB = SNR_dB-log(sps)/log(10)*10;  
%SNR_dB原始信噪比，log(sps)/log(10)*10将sps转换为分贝，SNR1_dB调整后的信噪比
%公式：SNR1_dB=SNR_dB?10?log10（sps）
P_sig = mean( abs(S1(:)).^2 );   
%S1(:)将S1展平为列向量，abs(S1(:)).^2计算信号每个采样点的功率，mean计算平均功率
P_noise = P_sig / (10^(SNR1_dB/10));
%10^(SNR1_dB/10)将信噪比从分贝转换为线性值，整体公式：根据信噪比计算噪声功率
if 0    
    % 如果是实值噪声:
    noise_matrix = sqrt(P_noise) * randn(1, length(Y));
    %randn(1,length(Y)):生成服从标准正态分布的随机数（方差为1）
    ...........................................................++++++++++
else
    % 如果是复值噪声:
    noise_matrix = sqrt(P_noise/2) * (randn(1, length(Y)) + 1i*randn(1, length(Y)));
    %sqrt(P_noise/2):调整噪声幅度，功率为P_noise;1i:虚数单位
end

Y_noise = Y + noise_matrix;
figure;plot(real(Y));
title('接受信号（含噪声）');

%% 按8192分段处理信号
segment_length=8192;
sign_length=length(Y_noise);
num_segments = floor(sign_length / segment_length)+1;

%% 
startflag=0;
for k = 1:num_segments
    % 提取当前段
    start_index = (k - 1) * segment_length + 1;  %当前段信号的起始点
    end_index = min(k * segment_length, sign_length);;    %当前段信号的终止点
    segment = Y_noise(start_index:end_index);   %截取指定段落的信号
    segment_Y = Y(start_index:end_index);   %截取的无噪声发射信号
    %% 能量检测法
    % 平方率检波
    segment_squared = abs(segment).^2;

    % 能量累积
    num_blocks = floor(length(segment_squared) / W);  % 能量块数量
    E = zeros(1, num_blocks);  % 能量块
    for lambda = 1:num_blocks
        start_idy = (lambda-1)*W + 1;
        end_idy = lambda*W;
        E(lambda) = sum(segment_squared(start_idy:end_idy));  % 能量累积
    end
    % 扩展能量累积信号
    ones1 = ones(W,1);
    E_spread = ones1*E;
    E_spread2 = E_spread(:).';

    % 调整 E_spread2 的长度与 Y 相同
    if length(E_spread2) < length(segment)
        E_spread2 = [E_spread2, zeros(1, length(segment) - length(E_spread2))];  % 补零
    else
        E_spread2 = E_spread2(1:length(segment));  % 截断
    end
    Esum = [Esum,  E_spread2];

    %% 降低信号Y与E

    % 降采样原始信号
    downsample_factor = 1;  % 降采样因子
    Y_downsampled = segment_Y(1:downsample_factor:end);  % 按间隔 downsample_factor 抽取采样点

    % 降采样能量累积信号
    E_spread2_downsampled = E_spread2(1:downsample_factor:end);



    %% 设置归一化门限
    Z_max = max(E);  % 能量序列的最大值
    Z_min = min(E);  % 能量序列的最小值
    gamma_norm = 0.5;  % 归一化门限比例
    gamma_facc = gamma_norm * (Z_max - Z_min) + Z_min;  % 固定门限
    %生成t时根据Y的长度动态调整
    t = (0:length(segment)-1) * T;  % 根据 Y 的长度生成时间向量

    %% 移动窗检测
    Wtemp = zeros(1,W_slip);
    for i =  1 : num_blocks-W_slip
        Wtemp =  E(i:i+W_slip-1);  %将能量赋值给移动窗
        Wtemp = Wtemp>gamma_facc;
        gamma_slip = sum(Wtemp);  %将窗能量求和
        %寻找起始点
        if  gamma_slip >  threadhold1 && startflag == 0%窗能量与阈值比较
            start_counter = start_counter + 1;
            toa_est = i ;  %定位起始点  
            startflag = 1;  %证明已经找到起始点，开始寻找结束点
            disp(['找到第',num2str(start_counter),'个起始点']);
            startpoint = [startpoint, i*W + (k-1)*segment_length];  %起始点数组
        end
        %寻找终止点
        if gamma_slip < W_slip-threadhold1 && startflag == 1
           finish_counter = finish_counter + 1;  % 终止点计数器加1
           toe_est = i;  %定位终止点
           startflag = 0;   %证明已经找到终止点，开始寻找起始点
           disp(['找到第 ', num2str(finish_counter), ' 个终止点']);
           finishpoint = [finishpoint, i*W + (k-1)*segment_length];  %终止点数组
        end
    end

end
    %% 绘制能量块结果
    figure;
    plot(Esum);  % 归一化能量累积信号
    hold on;
    plot(real(Y_noise));  % 归一化原始信号
    title('能量累积信号与原始信号 ');
    xlabel('采样点');
    ylabel('归一化幅度');
    legend('能量累积信号', '原始信号');
%% 结果
    figure;plot(real(Y_noise));  
    hold on;
    plot(startpoint, real(Y_noise(toa_est)), 'ro');
    plot(finishpoint,real(Y_noise(toa_est)), 'gs');
    title('信号起始终止点检测');
    xlabel('采样数');
    legend('接收信号', '估计的 TOA', '估计的 TOE');