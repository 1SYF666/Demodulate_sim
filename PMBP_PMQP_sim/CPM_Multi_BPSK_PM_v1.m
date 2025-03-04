%% ********* multi-carrier composite signals simulation ********* %%
%% ***** data:20250303 authoor:ShenYifu ****  %%
%{
    参考论文《复合调制信号盲分析技术研究》
    复合调制信号的识别与处理研究_2.1.3节多副载波复合调制信号
    增加了多载波的情况，查看外调制解调效果。
%}
%%
clc;clear;
close all;
%% 参数设置
K = 1e3;             % 单位 KHz
M = 1e6;             % 单位 MHz
Rb = 1*K;            % 码速率
fc = 1*K;            % 主载频
fcsubmax = 64*K;        % 副载频
fs = 5 * fcsubmax;      % 采样率
time = 5;            % 仿真时间单位秒
SNR = 20;

%% 多载波参数设置
Multinum = 3;
fcsub0 = fcsubmax/4;
Rbsub0 = 2*Rb;
fcsub1 = fcsubmax/2;
Rbsub1 = 2*Rb;
%% multi_BPSK_PM调制
fcsub = [fcsubmax,fcsub1,fcsub0];
Rbsub = [Rb,Rbsub1,Rbsub0];
N = fs * time;          % 总共采样点数
t = (0 : N-1)/fs;
s_BPSK =zeros(1,N);
for i = 1 : Multinum
    % PCM编码
    symbolnum = Rbsub(i)*time; % 码元个数
    pcm_bits = randi([0 1],1,symbolnum);
    pcm_symbols = 2 * pcm_bits - 1;   % 将0 --> -1; 1 --> 1
    sps = fs / Rbsub(i);       % 每符号采样数,取floor()防止非整数倍报错
    if rem(sps, 1) ~= 0
        disp("参数设置error,sps为小数");
        return;
    end
    rcos_pcm_s1 = repmat(pcm_symbols,sps,1);
    rcos_pcm_s2 = rcos_pcm_s1(:)';
    rcos_pcm_s = rcos_pcm_s2;
    % 多载波BPSK 调制
    s_BPSK = s_BPSK + rcos_pcm_s .* sin(2 * pi * fcsub(i) * t);
end
s_BPSK = s_BPSK/Multinum;  % 功率归一化

% PM 调制
Kp = pi/2;
s_PCM_BPSK_PM = cos(2 * pi * fc * t + Kp * s_BPSK);
s_PCM_BPSK_PM1 = sin(2 * pi * fc * t + Kp * s_BPSK);

send_signal = s_PCM_BPSK_PM + 1i * s_PCM_BPSK_PM1;
demodulation_signal = send_signal;
% figure;plot(fftshift(abs(fft(send_signal,65536))./max(abs(fft(send_signal,65536)))))
%% 信道
% 加噪
send_signal = awgn(send_signal, SNR, 'measured');

%% 测试

%%
send_signal = send_signal./max(abs(send_signal));
demodulation_signal = send_signal;
%% 参数估计
% 主载频估计
nfft = 16384*4;
fft_temp = abs(fft(send_signal,nfft));
fft_temp(1) = 0; % 去直流分量
[maxvalue,maxindex] = max(fft_temp);          % 选择全部

if  maxindex> nfft/2
    fc_est=(maxindex-nfft-1)/nfft*fs;% fc 为负值
else
    fc_est=(maxindex-1)/nfft*fs;  % fc 为正值
end
fprintf("主载频估计：%.3f\n",fc_est);

% 下变频
send_signal_temp = send_signal.*exp(-1i*2*pi*fc_est*(1:length(send_signal))/fs);
send_signal = send_signal_temp;

fcsub_est = fcsubmax;

%% 解调
% 主载波下变频
designla = demodulation_signal;
downsignal_temp = designla.*exp(-1i*2*pi*fc_est*(1:length(designla))/fs);
designla = downsignal_temp;

% PM解调
Signal_Channel_1 = designla;
Simulation_Length_1=length(Signal_Channel_1);

%参数置零
Signal_PLL_1=zeros(1,Simulation_Length_1);
NCO_Phase_1 =zeros(1,Simulation_Length_1);
Discriminator_Out_1=zeros(1,Simulation_Length_1);
Freq_Control_1=zeros(1,Simulation_Length_1);
PLL_Phase_Part_1=zeros(1,Simulation_Length_1);
PLL_Freq_Part_1=zeros(1,Simulation_Length_1);
I_PLL_1=zeros(1,Simulation_Length_1);
Q_PLL_1=zeros(1,Simulation_Length_1);

sigma = 0.707;
fs_nco = fs;
for i = 1:Simulation_Length_1
    coefficient_temp(i) = 0.01;
end
band_pm_bp = 2*fcsub_est;
BL=coefficient_temp*band_pm_bp;              
Wn=8*sigma*BL/(1+4*sigma^2);     T_nco=1/fs_nco;
K1(1:Simulation_Length_1)=(2*sigma*Wn*T_nco);
K2(1:Simulation_Length_1)=((T_nco*Wn).^2);

for i=2:Simulation_Length_1
    Signal_PLL_1(i)=Signal_Channel_1(i)*exp(-1i*(NCO_Phase_1(i-1)));
    I_PLL_1(i)=real(Signal_PLL_1(i));
    Q_PLL_1(i)=imag(Signal_PLL_1(i));
    Discriminator_Out_1(i) =(atan2( Q_PLL_1(i),I_PLL_1(i)));
    PLL_Phase_Part_1(i)=Discriminator_Out_1(i)*K1(i);
    PLL_Freq_Part_1(i)=Discriminator_Out_1(i)*K2(i)+PLL_Freq_Part_1(i-1);
    Freq_Control_1(i)=PLL_Phase_Part_1(i)+PLL_Freq_Part_1(i);
    NCO_Phase_1(i)=NCO_Phase_1(i-1)+Freq_Control_1(i)*2*pi;
end
figure;subplot(2,1,1);plot(fftshift(abs(fft(PLL_Phase_Part_1,65536)))); title("多载波频谱");
subplot(2,1,2);plot(fftshift(abs(fft(hilbert(PLL_Phase_Part_1),65536)))); title("多载波频谱");
PLL_Phase_Part_1 = hilbert(PLL_Phase_Part_1);     % 频谱上看，希尔伯特变换只保留了一半，使得后面下变频不用滤波
