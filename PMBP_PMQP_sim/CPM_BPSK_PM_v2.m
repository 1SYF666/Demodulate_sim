%% ********* CPM_BPSK_PM simulation ********* %%
%% ***** data:20250114 authoor:ShenYifu ****  %%
%% 参考论文《复合调制信号盲分析技术研究》
%{

%}
%%
clc;clear;
close all;
%% 参数设置
K = 1e3;             % 单位 KHz
M = 1e6;             % 单位 MHz
Rb = 1*K;            % 码速率
fc = 1*K;            % 主载频
fcsub = 64*K;        % 副载频
fs = 5 * fcsub;      % 采样率
time = 5;            % 仿真时间
symbolnum = Rb*time; % 码元个数
SNR = 20;

%% CPM_BPSK_PM调制
% PCM编码
pcm_bits = randi([0 1],1,symbolnum);
pcm_symbols = 2 * pcm_bits - 1;   % 将0 --> -1; 1 --> 1

sps = fs / Rb;                    % 每符号采样数
rcos_pcm_s1 = repmat(pcm_symbols,sps,1);
rcos_pcm_s2 = rcos_pcm_s1(:)';
rcos_pcm_s = rcos_pcm_s2;

% BPSK 调制
N = length(rcos_pcm_s);
t = (0 : N-1)/fs;
s_BPSK = rcos_pcm_s .* sin(2 * pi * fcsub * t);
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

%　副载波估计
M = 2;
s_suqre = send_signal.^M;
fft_temp = abs(fft(s_suqre,nfft));
fft_temp(1) = 0; % 去直流分量
[~,maxindex] = max(fft_temp);
if  maxindex> nfft/2
    fcsub_est=(maxindex-nfft-1)/nfft*fs/M;
else
    fcsub_est=(maxindex-1)/nfft*fs/M;
end
fcsub_est = abs(fcsub_est);
fprintf("副载频估计：%.3f\n",fcsub_est);

%　小波变换估计符号速率
cwtsignal = send_signal(1:65536);
[wt, freq, t,coi] = cwt(cwtsignal, 'amor', fs);

[idxvalue, idxClosest] = min(abs(freq - fcsub_est));
fprintf("最小相差：%.2f\n",idxvalue);
cwtfreal = wt(idxClosest, :,1);
cwtfrealabs = abs(cwtfreal);
% figure;plot(cwtfrealabs);title("小波系数图");

estsignal = cwtfrealabs;
meanest = mean(estsignal);
estsignalabs = abs(estsignal - meanest);
maxval1 = max(estsignalabs);
threadhold = 0.7*maxval1;
% figure;plot(estsignalabs);title("小波系数图");

% 最小间隔求符号速率
peakstemp = [];
flag = 1;
index = 1;
len = length(estsignalabs);
for i = 1 : len -1
    if flag  && estsignalabs(i)>threadhold && estsignalabs(i) > estsignalabs(i+1)
        flag = 0;
        peakstemp(index) = i;
        index = index +1;
    end

    if ~flag && estsignalabs(i) < estsignalabs(i+1)
        flag = 1;
    end
end

[mininterval,~] = min(abs(diff(peakstemp)));
Rb_est = fs/mininterval;
fprintf("符号速率估计：%.2f\n",Rb_est);


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
% figure;subplot(2,1,1);plot(fftshift(abs(fft(PLL_Phase_Part_1,65536)))); title("bpsk频谱");
% subplot(2,1,2);plot(fftshift(abs(fft(hilbert(PLL_Phase_Part_1),65536)))); title("bpsk频谱");
PLL_Phase_Part_1 = hilbert(PLL_Phase_Part_1);     % 频谱上看，希尔伯特变换只保留了一半，使得后面下变频不用滤波

% bpsk下变频
bpsksignal = PLL_Phase_Part_1;
downsignal_temp = bpsksignal.*exp(-1i*2*pi*fcsub_est*(1:length(bpsksignal))/fs);

% 重采样 为符号速的4倍
sps = 4;
fs_orignal = fs;
fs = sps*Rb_est;
[P,Q] = rat(fs/fs_orignal);
bpsksignal = resample(downsignal_temp, P, Q);

% bpsk载波同步
s3 = bpsksignal(1:end)/max(abs(bpsksignal)); % 输入
sum_number=length(s3);
Signal_Channel=zeros(1,sum_number);
Signal_PLL_D=zeros(1,sum_number);
NCO_Phase =zeros(1,sum_number);
Discriminator_Out=zeros(1,sum_number);Freq_Control=zeros(1,sum_number);
PLL_Phase_Part=zeros(1,sum_number);PLL_Freq_Part=zeros(1,sum_number);

Ko=1;Kd=1;Kp=Ko*Kd*1;
sigma=0.707;
T_nco=1/fs;
band_bpsk = Rb_est*2;          
BL = zeros(1,sum_number);Wn = zeros(1,sum_number);
C11 = zeros(1,sum_number);C22 = zeros(1,sum_number);
for i = 1:sum_number
    coefficient_temp=0.01;
    BL(i)=coefficient_temp*band_bpsk;
    Wn(i)=8*sigma*BL(i)/(1+4*sigma^2);     
    C11(i)=(2*sigma*Wn(i)*T_nco)/(Kp);      
    C22(i)=((T_nco*Wn(i))^2)/(Kp);
end

Vi_PLL=s3;
Vi_PLL_real=real(Vi_PLL);  Vi_PLL_imag=imag(Vi_PLL);
for i=2:sum_number
    Signal_Channel(i)=Vi_PLL_real(i)+1i*Vi_PLL_imag(i);
    Signal_PLL_D(i)=Signal_Channel(i)*exp(-1i*(mod(NCO_Phase(i-1),2*pi)));
    Vo_PLL_real(i)=real(Signal_PLL_D(i));   
    Vo_PLL_imag(i)=imag(Signal_PLL_D(i));
    Discriminator_Out(i)=sign(Vo_PLL_real(i))*Vo_PLL_imag(i);
    PLL_Phase_Part(i)=Discriminator_Out(i)*C11(i);
    PLL_Freq_Part(i)=Discriminator_Out(i)*C22(i)+PLL_Freq_Part(i-1);
    Freq_Control(i)=PLL_Phase_Part(i)+PLL_Freq_Part(i);
    NCO_Phase(i)=NCO_Phase(i-1)+Freq_Control(i)*2*pi;
end
%
% figure;plot(PLL_Freq_Part*fs);title("载波同步频率跟踪曲线");
% figure;plot(real(exp(-1i*(mod(NCO_Phase,2*pi)))));title("载波输出曲线");
% figure;plot(Vo_PLL_real);title("载波同步输出:real部");
% figure;plot(Vo_PLL_imag);title("载波同步输出:imag部");
% figure;scatter(Vo_PLL_real(5000:end),Vo_PLL_imag(5000:end));title("载波同步输出星座图");

% bpsk符号同步
s_bpsk = Vo_PLL_real + 1i*Vo_PLL_imag;
aI=Vo_PLL_real(1:end);
bQ=Vo_PLL_imag(1:end);
N=floor(length(aI)/4); %符号数  floor向负无穷取整
Ns=4*N;  %总的采样点数
w=[0.5,zeros(1,N-1)];  %环路滤波器输出寄存器，初值设为0.5
n=[0.9 zeros(1,Ns-1)]; %NCO寄存器，初值设为0.9  可调
n_temp=[n(1),zeros(1,Ns-1)];
u=[0.6,zeros(1,2*N-1)];%NCO输出的定时分数间隔寄存器，初值设为0.6
yI=zeros(1,2*N);       %I路内插后的输出数据
yQ=zeros(1,2*N);       %Q路内插后的输出数据
time_error=zeros(1,N); %Gardner提取的时钟误差寄存器
ik=time_error;
qk=time_error;
i=1;    %用来表示Ts的时间序号,指示n,n_temp,nco,
kk=1;   %用来表示Ti时间序号,指示u,yI,yQ
ms=1;   %用来指示T的时间序号,用来指示a,b以及w
strobe=zeros(1,Ns);
BL1=Rb_est*0.01;                %可调
decimator=1;
Wn1=8*sigma*BL1/(1+4*sigma^2);    %环路自由震荡角频率
T_nco1=1/(fs)*decimator;   %压控振荡器NCO频率更新周期fs_demo/K_sampdemo
c1=(2*sigma*Wn1*T_nco1);      %环路滤波器系数c1
c2=((T_nco1*Wn1)^2);          %环路滤波器系数c2
ns=length(aI)-2;
while(i<ns)
    n_temp(i+1)=n(i)-w(ms);

    if(n_temp(i+1)>0)
        n(i+1)=n_temp(i+1);
    else
        n(i+1)=(n_temp(i+1)-ceil(n_temp(i+1))+1);
        %内插滤波器模块
        FI1=0*aI(i+2)-1/2*aI(i+1)+1/2*aI(i)-1/6*aI(i-1);
        FI2=1/6*aI(i+2)+1/2*aI(i+1)-1*aI(i)+1/2*aI(i-1);
        FI3=-1/6*aI(i+2)+1*aI(i+1)-1/2*aI(i)-1/3*aI(i-1);
        FI4=0*aI(i+2)+0*aI(i+1)+1*aI(i)+0*aI(i-1);
        yI(kk)=(((FI1*u(kk)+FI2)*u(kk)+FI3)*u(kk)+FI4);

        FQ1=0*bQ(i+2)-1/2*bQ(i+1)+1/2*bQ(i)-1/6*bQ(i-1);
        FQ2=1/6*bQ(i+2)+1/2*bQ(i+1)-1*bQ(i)+1/2*bQ(i-1);
        FQ3=-1/6*bQ(i+2)+1*bQ(i+1)-1/2*bQ(i)-1/3*bQ(i-1);
        FQ4=0*bQ(i+2)+0*bQ(i+1)+1*bQ(i)+0*bQ(i-1);
        yQ(kk)=(((FQ1*u(kk)+FQ2)*u(kk)+FQ3)*u(kk)+FQ4);

        strobe(kk)=mod(kk,2);
        %时钟误差提取模块，采用的是Gardner算法
        if(strobe(kk)==0)
            %取出插值数据
            ik(ms)=yI(kk);
            qk(ms)=yQ(kk);
            %每个数据符号计算一次时钟误差
            if(kk>2)
                Ia=(yI(kk)+yI(kk-2))/2;
                Qa=(yQ(kk)+yQ(kk-2))/2;
                time_error(ms)=(yI(kk-1)-Ia)*(yI(kk)-yI(kk-2))+(yQ(kk-1)-Qa)*(yQ(kk)-yQ(kk-2));
            else
                time_error(ms)=(yI(kk-1)*yI(kk)+yQ(kk-1)*yQ(kk));
            end
            %环路滤波器,每个数据符号计算一次 环路滤波器输出
            if(ms>1)
                w(ms+1)=w(ms)+c1*(time_error(ms)-time_error(ms-1))+c2*time_error(ms-1);
            else
                w(ms+1)=w(ms)+c1*time_error(ms)+c2*time_error(ms);
            end
            ms=ms+1;
        end
        kk=kk+1;
        u(kk)=n(i)/w(ms);
    end

    i=i+1;
end
I_PLL_D1=ik(1:end);
Q_PLL_D1=qk(1:end);

% figure;plot(u);title("小数时间间隔u,符号数的2倍");
% figure;plot(n);title("NCO寄存器内容n,符号数的4倍");
% figure;plot(w(1:end-1));title("经过环路滤波器得到的定时控制字w，符号数");
% figure;plot(time_error);title("定时误差估计值，符号数");
%% 绘制星座图
starmap = (I_PLL_D1+1i*Q_PLL_D1).*exp(-1i*0.25*pi);
figure;scatter(real(starmap(100:end)),imag(starmap(100:end)));title("码元同步输出星座图");