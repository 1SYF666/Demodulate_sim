%% ********* DDC_F simulation ********* %%
%% ***** data:20250126 authoor:ShenYifu ****  %%
%%
%{
    参考论文：
    1.《基于FPGA的FM调制解调器的实现_曹沅》
    2.《基于软件无线电的数字调制解调技术的研究与实现_王松涛》

%}
%%
clear; clc;
close all;
%% 读文件
% str = 'D:\MODIFI_CODE\协议信号研究\A1\DDC_F225MHz_T2022-07-17_17-13-50-972_00400000sps_CH2.std';    % 2025-01-27
str = 'D:\MODIFI_CODE\协议信号研究\A1\DDC_F225MHz_T2022-07-17_17-14-17-972_00400000sps_CH2.std';    % 2025-01-28

fileID = fopen(str,'r');
data = fread(fileID,'int16');
fclose(fileID);
data = data.';
s = zeros(1,2*length(data));
s = data(1:2:end-1) + 1i * data(2:2:end);

%% 截取信号部分
if 1
    % 第一部分
    % signalfirstsegmentaion = s(400000:449142);   % 长度参考cooledit pro软件设置
    % signalfirstsegmentaion = s(4151116:4184266);   % 2025-01-28
else
    % 第二部分
    % signalfirstsegmentaion = s(510000:584435);  % 2025-01-27
    signalfirstsegmentaion = s(4276617:4328711);  % 2025-01-28
end
%
signalfirstsegmentaion = signalfirstsegmentaion./max(abs(signalfirstsegmentaion));

%% 参数设置
fs = 400e3;


%% 绘图
if 0
    figure;title('时域和频域波形');subplot(2,1,1);plot(real(s));subplot(2,1,2);plot(imag(s));
    % 绘制时间频率图（频谱图）
    figure;
    spectrogram(s, 65536, 16384, 65536, fs, 'yaxis');
    title('Time-Frequency Representation (Spectrogram)');
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
end

%% 信号分析部分

if 0
    % 绘制功率谱
    % 使用 periodogram 计算功率谱
    [pxx, f] = periodogram(signalfirstsegmentaion, hamming(length(signalfirstsegmentaion)), length(signalfirstsegmentaion), fs);

    % 绘制功率谱
    figure;
    plot(f, 10*log10(pxx));
    grid on;
    title('Power Spectral Density using Periodogram');
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');

    % 绘制包络谱
    s_envelope = real(signalfirstsegmentaion).^2 + imag(signalfirstsegmentaion).^2 ;
    envelope = abs(hilbert(signalfirstsegmentaion));  % Hilbert 变换获取信号包络
    figure;plot(envelope); title("包络谱示意图");
    figure;plot(s_envelope); title("包络谱示意图");
    figure;plot(fftshift(abs(fft(s_envelope,65536)))); title("包络谱示意图");
end

% 分析信号感觉像FM调制
% 载频估计
send_signal = signalfirstsegmentaion(1:end);
nfft = 16384*4*4;
fft_temp = abs(fft(send_signal,nfft));
fft_temp(1) = 0; % 去直流分量
[maxvalue,maxindex] = max(fft_temp);          % 选择全部

if  maxindex> nfft/2                                 % 后面运行此模块代码发现负频率估计失效，
    % fc 为负值                                      % 故又添加了此处的if判断语句
    fc_est=(maxindex-nfft-1)/nfft*fs;
else
    % fc 为正值
    fc_est=(maxindex-1)/nfft*fs;
end
fprintf("载频估计：%.3f\n",fc_est);

% FM解调
if 1
    s_MSK_I = real(send_signal);
    s_MSK_Q = imag(send_signal);
    for k = 13 : 13
        for i = k+1:1:length(s_MSK_I)
            Pdot(i-1) = s_MSK_I(i-k)*s_MSK_I(i) + s_MSK_Q(i-k)*s_MSK_Q(i);
            Pcross(i-1) = s_MSK_I(i-k)*s_MSK_Q(i) - s_MSK_Q(i-k)*s_MSK_I(i);
        end
        theta_w = atan2(Pcross,Pdot);
        theta_w = theta_w - mean(theta_w);
        figure;subplot(3,1,1);plot(theta_w);title('鉴相器输出');subplot(3,1,2);plot(fftshift(abs(fft(theta_w))));title('鉴相器输出频谱图');
        subplot(3,1,3);scatter(theta_w,theta_w);title('鉴相器输出星座图');
        % figure;cwt(theta_w, 'amor', fs);

        if 0
            % 滤波器
            hai = 3.3;                             %海明窗窗过度带宽系数
            wp = 0.16*pi;                          %通带截止频率 wp < fc/(2*fs)
            ws = 0.26*pi;                           %阻带起始频率 ws > fc/(2*fs)
            wdlta = ws-wp;
            N_lp = ceil(2*pi*hai/wdlta);           %求滤波器阶数N_lp
            Wc = (wp+ws)/2;
            b = fir1(N_lp-1,Wc/pi,hamming(N_lp));
            theta_w_filter1 = filter(b,1,theta_w);
            figure;subplot(2,1,1);plot(theta_w_filter1);title('滤波后鉴相器输出');subplot(2,1,2);plot(fftshift(abs(fft(theta_w_filter1))));
            theta_w = theta_w_filter1;
        end

        % 符号速率估计
        estrbsignal = theta_w( 1:end );

        nfft  = 16384*4;
        M =  2;
        s_envalop = estrbsignal.*conj(estrbsignal);
        s_temp = diff(s_envalop);
        fft_temp = abs(fft(s_temp,nfft));
        fft_temp(1:20) = 0;     % 消除干扰
        half_fft_temp = fft_temp(1:nfft/2);
        figure;plot(half_fft_temp);
        [maxvalue,maxindex] = max(half_fft_temp);
        Rb_est = ( (maxindex-1)*fs/nfft);
        fprintf("码速估计：%.3f\n",Rb_est);
    end
else
    s_re_filter1 = real(send_signal(1:end));
    s_im_filter1 = imag(send_signal(1:end));

    % 延迟相乘
    minlength = min(length(s_re_filter1),length(s_im_filter1));

    for k = 13 : 13
        for i = 1 : minlength-k-1
            theta_w(i) = s_re_filter1(i)*s_im_filter1(i+k)-s_im_filter1(i)*s_re_filter1(i+k);
        end
        figure;subplot(3,1,1);plot(theta_w);title('鉴相器输出');subplot(3,1,2);plot(fftshift(abs(fft(theta_w))));
        subplot(3,1,3);scatter(theta_w,theta_w);
    end
end

% 同步
if 0
    % 搜索起始点
    diffsignal = abs(diff(theta_w));
    % figure;plot(diffsignal);
    segementaiondiff = diffsignal(1: 100);
    [valuemax,indexmax]= max(segementaiondiff);
    startpoint = indexmax + k;

    % 降采样
    deout = theta_w(startpoint:k:end-1);
    figure;subplot(3,1,1);plot(deout);title('鉴相器降采样输出');subplot(3,1,2);plot(fftshift(abs(fft(deout))));
    subplot(3,1,3);scatter(deout,deout);
else
    % 重采样 为符号速的4倍
    Rb_est = fs/(k*2);
    sps = 4;
    fs_orignal = fs;
    fs = sps*Rb_est;
    [P,Q] = rat(fs/fs_orignal);
    theta_wsignal = resample(theta_w, P, Q);
    figure;subplot(3,1,1);plot(theta_wsignal);title('鉴相器降采样输出');subplot(3,1,2);plot(fftshift(abs(fft(theta_wsignal))));
    subplot(3,1,3);scatter(theta_wsignal,theta_wsignal);

    theta_wsignal = theta_wsignal./max(abs(theta_wsignal));
    Vo_PLL_real = theta_wsignal;
    Vo_PLL_imag = theta_wsignal;

    % 符号同步
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
    decimator=1;sigma = 0.707;
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

    figure;subplot(4,1,1);plot(u);title("小数时间间隔u,符号数的2倍");
    subplot(5,1,2);plot(n);title("NCO寄存器内容n,符号数的4倍");
    subplot(5,1,3);plot(w(1:end-1));title("经过环路滤波器得到的定时控制字w，符号数");
    subplot(5,1,4);plot(time_error);title("定时误差估计值，符号数");
    subplot(5,1,5);scatter(I_PLL_D1(200:end-1),Q_PLL_D1(200:end-1));title("码元输出星座图");

end


















