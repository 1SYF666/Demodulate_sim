%% ************* wcdma sigle user sim ************ %%
%% *******  date:20241129 author:ShenYifu ******* %%
%% *********************************************** %%
%{
    note:
        1. 对此前代码进行更新，使其更加规范
        2. 更新解扩代码-迟早门算法
        3. 对10ms数据进行处理，相当于一个突发数据
           Rb_I * fream_t * spread_factor_I = length(Gold)
           1920000 * 0.01 * 2 = 38400
%}
%{
    注意：解扰解扩分开进入迟早门
%}
%%
clc;clear;
close all;

%% 参数设置
Rb = 384000;    
pllloopfalg = 0;                    % 鉴相器标志位
if Rb == 600
    spread_factor_I = 256;
elseif Rb == 2400
    spread_factor_I = 256;
elseif Rb == 32000
    spread_factor_I = 32;
elseif Rb == 64000
    spread_factor_I = 16;
    pllloopfalg = 1;
elseif Rb == 128000
   spread_factor_I = 8;
   pllloopfalg = 1;
elseif Rb == 384000
    spread_factor_I=4;
    % spread_factor_I=2;
    pllloopfalg = 1;
end
spread_factor_Q = 256;
Rb_scarmble = 3.84e6;
sf_num = 2;
Rb_I = Rb_scarmble/spread_factor_I;           
Rb_Q = Rb_scarmble/spread_factor_Q;           
sps = 5;
fs = Rb_I*spread_factor_I*sps;            
fream_t = 0.01;           % 10ms                               
SNR = 5;                                
%% 信号生成
[signal,Gold,ovsf_i,ovsf_q,...
despread_i,despread_q,data_code_A,source_A]...
=single_wcdma_generate(fream_t,Rb,Rb_I,spread_factor_I,sf_num,SNR,sps);


%% 测试载频影响 
% Vi_est_down=signal.*exp(-2*pi*(100)*( (1:1:length(signal))/(sps*3.84e6) )*sqrt(-1));
Vi_est_down = signal;

%% 捕获参数初始化
acquisition_len = 8192;
chip_len = acquisition_len * sps;
gold_chip =  Gold(1:acquisition_len);
acquisition_sequence = repmat(gold_chip,sps,1);
acquisition_sequence = acquisition_sequence(:).';
acquisition_fft = conj (fft(acquisition_sequence,chip_len));
Vi_acquisition = Vi_est_down;

%% 捕获
chip_num = floor(length(Vi_est_down)/chip_len);
flag = 0;
for i = 1 : chip_num
    acquisition_temp = Vi_acquisition((1+chip_len*(i-1)):(chip_len*i));
    autocorrelation = abs( ifft(fft(acquisition_temp).*acquisition_fft) );
    [maxtemp(1,1),maxtemp(2,1)] = max(autocorrelation);
    if maxtemp(1,1)>3000&&flag == 0
        delay_est = maxtemp(2,1)+(i-1)*chip_len;
        % plot(autocorrelation);
        flag = 1;
        break;
    end
end

s_capture = Vi_acquisition(delay_est+2:end);
fprintf("延迟估计：\t\t\t%d\n",delay_est);

%% 解扰初始化操作
Vi_track = s_capture;

% 扰码上采样
gold_sps = repmat(Gold,sps,1);
gold_sps_conj = gold_sps(:)';     % 扰码的共轭，用于解扩

descramblelen = spread_factor_Q*sps;     %　一次解扰输出长度                        
descramble_s = Vi_track;
estPN_sps_acq = gold_sps_conj;
num_symbol = floor(length(descramble_s)/(descramblelen));

% 环路初始化
looprb = Rb_Q;
loopinitial;
%% 解扰操作
for i = 1:num_symbol
    receivetemp = descramble_s( 1 + descramblelen * (i-1) : descramblelen * i );

    Early = circshift(estPN_sps_acq,-floor(sps/2));
    Early = Early(1 + descramblelen * (i-1) : descramblelen * i);
    Prompt = circshift(estPN_sps_acq, 0);
    Prompt = Prompt(1 + descramblelen * (i-1) : descramblelen * i);
    Late = circshift(estPN_sps_acq,floor(sps/2));
    Late = Late(1 + descramblelen * (i-1) : descramblelen * i);
    
    % 超前
    i_E = receivetemp.*Early;
    I_E = sum(i_E); 
    E(i) = abs(sqrt(I_E.^2));
    % 即时
    i_P = receivetemp.*Prompt;
    I_P = sum(i_P);    
    P(i) = abs(sqrt(I_P.^2));
    output_tra(1 + descramblelen * (i-1) : descramblelen * i) = i_P;
    % 滞后
    i_L = receivetemp.*Late; 
    I_L = sum(i_L);
    L(i) = abs(sqrt(I_L.^2));
    % 环路--测试对比使用
    j = i+1;
    loop_error_In(j) = 1/2*(E(i)^2 - L(i)^2)/(E(i)^2 + L(i)^2);
    PLL_Phase_Part(j)=loop_error_In(j)*C11;
    PLL_Freq_Part(j)=loop_error_In(j)*C22+PLL_Freq_Part(j-1);
    loop_error_Out(j)=PLL_Phase_Part(j)+PLL_Freq_Part(j);
    temp2_diff = loop_error_Out(j)-loop_error_Out(j-1);

    %    if i>1000
    % 比较EPL大小
    if L(i)>P(i)&&E(i)<P(i)
        % 信号滞后PN码一个点
        estPN_sps_acq_temp = circshift(estPN_sps_acq,floor(sps/4));
        estPN_sps_acq = estPN_sps_acq_temp;
    end
    if L(i)<P(i)&&E(i)>P(i)
        % 信号超前PN码一个点
        estPN_sps_acq_temp = circshift(estPN_sps_acq,-floor(sps/4));
        estPN_sps_acq = estPN_sps_acq_temp;
    end
    % 还需考虑失锁的情况
end
descramblesignal = output_tra;
%% 解扩初始化操作
despreadlen = spread_factor_Q*sps;     %　一次解扩输出长度  

% i路扩频码上采样
spreadmultiple = spread_factor_Q/spread_factor_I ;
ovsf_i_temp = repmat(ovsf_i,1,spreadmultiple);
ovsf_i_sps = repmat(ovsf_i_temp,sps,1);
ovsf_i_sps = ovsf_i_sps(:).';

% q路扩频码上采样
ovsf_q_sps = repmat(ovsf_q,sps,1);
ovsf_q_sps = ovsf_q_sps(:).';
pn_sps_i = ovsf_i_sps;
pn_sps_q = ovsf_q_sps;

r_dec6 = descramblesignal;
num_symbol = floor(length(descramblesignal)/(despreadlen));
% 环路初始化
looprb = Rb_Q;
loopinitial;
%% 解扩操作
for i = 1:num_symbol
    rr_dec6 = r_dec6( 1 + despreadlen * (i-1) : despreadlen * i );
    rr_dec6_I = real(rr_dec6);
    rr_dec6_Q = imag(rr_dec6);

    Earlyi = circshift(pn_sps_i,-floor(sps/2));
    Earlyq = circshift(pn_sps_q,-floor(sps/2));
    Prompti = circshift(pn_sps_i, 0);
    Promptq = circshift(pn_sps_q, 0);
    Latei = circshift(pn_sps_i,floor(sps/2));
    Lateq = circshift(pn_sps_q,floor(sps/2));

    % 超前
    i_E = rr_dec6_I.*Earlyi;q_E = rr_dec6_Q.*Earlyq;
    I_E = sum(i_E); Q_E = sum(q_E);
    E(i) = abs(sqrt(I_E.^2 + Q_E.^2));
    % 即时
    i_P = rr_dec6_I.*Prompti; q_P = rr_dec6_Q.*Promptq;
    I_P = sum(i_P);Q_P = sum(q_P);
    P(i) = abs(sqrt(I_P.^2+Q_P.^2));
    output_tra(1 + despreadlen * (i-1) : despreadlen * i) = i_P + 1i * q_P;
    % 滞后
    i_L = rr_dec6_I.*Latei; q_L = rr_dec6_Q.*Lateq;
    I_L = sum(i_L); Q_L = sum(q_L);
    L(i) = abs(sqrt(I_L.^2+Q_L.^2));
    % 环路
    j = i+1;
    loop_error_In(j) = 1/2*(E(i)^2 - L(i)^2)/(E(i)^2 + L(i)^2);
    PLL_Phase_Part(j)=loop_error_In(j)*C11;
    PLL_Freq_Part(j)=loop_error_In(j)*C22+PLL_Freq_Part(j-1);
    loop_error_Out(j)=PLL_Phase_Part(j)+PLL_Freq_Part(j);
    temp2_diff = loop_error_Out(j)-loop_error_Out(j-1);

    % 比较EPL大小
    if L(i)>P(i)&&E(i)<P(i)
        % 信号滞后PN码一个点
        pn_sps_i_temp = circshift(pn_sps_i,floor(sps/4));
        pn_sps_i = pn_sps_i_temp;
        pn_sps_q_temp = circshift(pn_sps_q,floor(sps/4));
        pn_sps_q = pn_sps_q_temp;

    end
    if L(i)<P(i)&&E(i)>P(i)
        % 信号超前PN码一个点
        pn_sps_i_temp = circshift(pn_sps_i,-floor(sps/4));
        pn_sps_i = pn_sps_i_temp;
        pn_sps_q_temp = circshift(pn_sps_q,floor(sps/4));
        pn_sps_q = pn_sps_q_temp;
    end
    
end
despreadsignal = output_tra;
figure ;plot(real(despreadsignal));title("I路");
figure;plot(imag(despreadsignal));title("Q路");

%% 载波同步初始化
Vo_des = despreadsignal;
sum_length = spread_factor_I*sps;
sum_number = floor(length(Vo_des)/(sum_length));
for i = 1 : sum_number
    Vo_sum(i) = sum(Vo_des((i-1)*sum_length+1:i*sum_length))/sum_length ; 
end
looprb = fs/sum_length;
if Rb == 384000
    loopcoef = 0.02;
elseif Rb == 128000
    loopcoef = 0.1;
elseif Rb == 64000
    loopcoef = 0.01;
else
    loopcoef = 0.001;
end

pll_loop_initial;  

if pllloopfalg == 1
    Vi_PLL_real = real(Vo_sum);  Vi_PLL_imag = real(Vo_sum);
else
    Vi_PLL_real = real(Vo_sum);  Vi_PLL_imag = imag(Vo_sum);
end
%% 载波同步
for i=2:Simulation_Length
    Signal_Channel(i)=Vi_PLL_real(i)+1i*Vi_PLL_imag(i);
    Signal_PLL_D(i)=Signal_Channel(i)*exp(-1i*(mod(NCO_Phase(i-1),2*pi)));    
    Vo_PLL_real(i)=real(Signal_PLL_D(i));  
    Vo_PLL_imag(i)=imag(Signal_PLL_D(i));
    if pllloopfalg == 1
        Discriminator_Out(i)=sign(Vo_PLL_real(i))*Vo_PLL_imag(i);
    else
        Discriminator_Out(i)=sign(Vo_PLL_real(i))*Vo_PLL_imag(i)-sign(Vo_PLL_imag(i))*Vo_PLL_real(i);
    end
    PLL_Phase_Part(i)=Discriminator_Out(i)*C11;
    PLL_Freq_Part(i)=Discriminator_Out(i)*C22+PLL_Freq_Part(i-1);
    Freq_Control(i)=PLL_Phase_Part(i)+PLL_Freq_Part(i);
    NCO_Phase(i)=NCO_Phase(i-1)+Freq_Control(i)*2*pi;
    
end
figure;plot(PLL_Freq_Part);title("载波同步频率跟踪曲线");

figure;plot(Vo_PLL_real);hold on;plot(1.2*(despread_i*2-1));title("解扩后与加扩前I路数据对比");

%% 解扩后误码率计算
reals = sign(Vo_PLL_real);
statrsymbol = 200;

% printf("")



