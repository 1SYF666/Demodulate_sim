%% ********* music simulation ********* %%
%% ***** data:20250311 authoor:ShenYifu ****  %%
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
SNR = 20;

nums = 1; % 信号源数
numy = 3; % 接收阵元数
d = 0.1;  % 单位米
theta = 20/180*pi; % 角度
c = 3*10^8;% 光速

%% BPSK信号生成
% PCM编码
pcm_bits = randi([0 1],1,symbolnum);
pcm_symbols = 2 * pcm_bits - 1;   % 将0 --> -1; 1 --> 1
% BPSK 调制
sps = fs / Rb;                    % 每符号采样数
rcos_pcm_s1 = repmat(pcm_symbols,sps,1);
rcos_pcm_s1(2:sps,:) = 0;
rcos_pcm_s = rcos_pcm_s1(:)';
rolloff_factor = 0.7;       % 滚降因子
rcos_fir = rcosdesign(rolloff_factor,2*sps,sps); % 默认是根升余弦滤波器,'sqrt'
rcos_ds_real = filter(rcos_fir,1,rcos_pcm_s);

%% 单源信号模型
S1 = rcos_ds_real;
fa = exp(-j*2*pi*fc*d*sin(theta)/c);
exponent = linspace(0,numy-1,numy)';
A = fa.^exponent;

Y = A * S1;

%% 信道

%% 接收端

% 空间协方差矩阵：
Ry = Y*Y'/length(Y); 

% 特征值分解
[U, Lambda] = eig(Ry); % Lambda 为对角阵，U 的列为特征向量

% 噪声子空间
[sortvalue,beforeindex] = sort(diag(Lambda),'descend');

U1 = U(:,beforeindex);
Un = U1(:,numy-nums:end);

% 构造估计阵列流形矩阵
theta_est = [0:0.5:179]/180*pi;
fa_est = exp(-j*2*pi*fc*d*sin(theta_est)/c);
A_est = fa_est.^exponent;

% 计算music伪谱
for k =  1: length(theta_est)
    Pmusic(k) = 1/(A_est(:,k)'*Un*Un'*A_est(:,k));
end
figure;plot(theta_est/pi*180,abs(Pmusic));title('传统方法music谱');





