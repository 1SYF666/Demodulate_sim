%% ********* music simulation ********* %%
%% ***** data:20250311 authoor:ShenYifu ****  %%
%{
参考论文：
        《基于特征值聚类的MUSIC算法》  
%}
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
SNR_dB = 30;        % (单位 dB);
numy = 10; % 接收阵元数
d = 0.1;  % 单位米
theta = [20/180*pi 40/180*pi 50/180*pi 60/180*pi]; % 角度
fc = 1500*M;
nums = length(theta); % 信号源数

c = 3*10^8;% 光速

%% BPSK信号生成
% PCM编码
pcm_bits = randi([0 1],nums,symbolnum);
pcm_symbols = 2 * pcm_bits - 1;   % 将0 --> -1; 1 --> 1

% BPSK 调制
sps = fs / Rb;                    % 每符号采样数
rolloff_factor = 0.7;       % 滚降因子
rcos_fir = rcosdesign(rolloff_factor,2*sps,sps); % 默认是根升余弦滤波器,'sqrt'
for i = 1:nums
    pcm_symbols_temp = pcm_symbols(i,:);
    rcos_pcm_s1 = repmat(pcm_symbols_temp,sps,1);
    rcos_pcm_s1(2:sps,:) = 0;
    rcos_pcm_s = rcos_pcm_s1(:)';
    rcos_ds_real(i,:) = filter(rcos_fir,1,rcos_pcm_s);
end

%% 多源信号模型
% S_coefficience = ones(nums,1);
% S1 = rcos_ds_real;
% S_multiple = S_coefficience*S1;
S_multiple = rcos_ds_real;
fa = exp(-j*2*pi.*fc*d.*sin(theta)/c);
exponent = linspace(0,numy-1,numy)';
A = fa.^exponent;
Y = A * S_multiple;
%% 信道
SNR1_dB = SNR_dB-log(sps)/log(10)*10;
P_sig = mean( abs(Y(:)).^2 ); 
P_noise = P_sig / (10^(SNR1_dB/10));
if 0
    % 如果是实值噪声:
    noise_std = sqrt(P_noise);
    noise_matrix = noise_std * randn(numy, length(Y));
else
    % 如果是复值噪声:
    noise_matrix = sqrt(P_noise/2) * (randn(numy, length(Y)) + 1i*randn(numy, length(Y)));
end

Y_noise = Y + noise_matrix;
%% 接收端
% 空间协方差矩阵：
Ry = Y_noise*Y_noise'/length(Y_noise); 

% 特征值分解
[U, Lambda] = eig(Ry); % Lambda 为对角阵，U 的列为特征向量

% 估计信源个数 -- 阈值可调节参考mydescan_threshold_v2.m
eps    = 0.3;    % 邻域半径
minPts = 5;      % 最小邻域点数
datadiag = diag(Lambda);
disp(['共有：',num2str(length(datadiag)),'个特征值']);
disp(['特征值：',num2str(datadiag')]);
labels = myDBSCAN(datadiag, eps, minPts);
nums_est = sum(labels==0);
disp(['估计信源个数:',num2str(nums_est)]);
% 噪声子空间
[sortvalue,beforeindex] = sort(diag(Lambda),'descend');
U1 = U(:,beforeindex);
Un = U1(:,numy-nums_est:end);

% 构造估计阵列流形矩阵
theta_est = [0:0.5:90]/180*pi;
fcest = fc;
fa_est = exp(-j*2*pi*fcest*d*sin(theta_est)/c);
A_est = fa_est.^exponent;

% 计算music伪谱
for k =  1: length(theta_est)
    Pmusic(k) = 1/(A_est(:,k)'*Un*Un'*A_est(:,k));
end
figure;plot(theta_est/pi*180,10*log10(abs(Pmusic)));title('传统方法music谱');
xlabel('角度'); ylabel('伪谱幅度');

% 估计角度
clear sortvalue beforeindex;
[sortvalue,beforeindex]= sort(abs(Pmusic));

theta_index = beforeindex(end-nums_est+1:end);
doa_theta = theta_est(theta_index)/pi*180;
disp(['估计角度为：',num2str(doa_theta)]);
disp(['估计偏差为：',num2str(doa_theta-theta/pi*180),'°']);
