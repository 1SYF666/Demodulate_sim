%{
    输入参数：
            fream_t--帧时长；fs--采样率；Rb--信息速率；Rb_I--I路信道速率；
            spread_factor_i--I路扩频码长；sf_num--扩频码编号；SNR--信噪比；
            Q路信息为控制信道，Q路信息、符号速率15K、以及扩频码长为256都是固定的。
    输出参数：
            send_signal--发射端输出信号；gold--扰码；ovsf_codes_i--I路扩频码；
            ovsf_codes_q--Q路扩频码；channel_data_i--信道编码后i路数据；
            channel_data_q--信道编码后q路数据；data_code -- CRC后信息数据；
            data_i -- 信息比特
    I路信息速率：600  2400 32000 64000 128000 384000 
%}
function [send_signal,gold,ovsf_codes_i,ovsf_codes_q,...
    channel_data_i,channel_data_q,data_code,data_i]...
    = single_wcdma_generate(fream_t,Rb,Rb_I,spread_factor_i,sf_num,snr,sps)  
%% 参数设置  ——————  需要加信道编码，如果我有时间的话
j = sqrt(-1);
% I路参数设置  
if Rb==600
    padding = 28+16;                
    trellis = poly2trellis(9, [771 663 557]);                   %（3，1，9）卷积编码
elseif Rb==2400
    padding = 10+16;                
    trellis = poly2trellis(9, [771 663 557]);                   %（3，1，9）卷积编码
elseif Rb==32000
    padding=64+16;
    trellis = poly2trellis(9, [771 663 557]);                   %（3，1，9）卷积编码
elseif Rb==64000
    padding=144+16;
    trellis = poly2trellis(9, [771 663 557]);                   %（3，1，9）卷积编码
elseif Rb==128000
   padding=304+16;
   trellis = poly2trellis(9, [771 663 557]);                    %（3，1，9）卷积编码
elseif Rb==384000
    if spread_factor_i==4
        trellis = poly2trellis(9, [561 753]);                     %（2，1，9）卷积编码
        padding=944+16;
    elseif spread_factor_i==2
        trellis = poly2trellis(9, [561 753]);
        padding=5760;
    end
end
sim_base_i = fream_t*Rb;                                        % 信息比特个数
sim_channel_i = fream_t*Rb_I;                                   % 信道比特个数
ovsf_codes_i = ovsf_generator(spread_factor_i,sf_num);          % I路扩频码

% Q路参数设置
Rb_Q = 15000;                                                   % Q路信道速率
spread_factor_q = 256;              
sim_channel_q = fream_t*Rb_Q;
ovsf_codes_q = ovsf_generator(spread_factor_q,1);               % Q路扩频码

%% 信息比特
data_i = randi([0 1],1,sim_base_i);                             % in 10 ms, 15 slots

%% 信道编码 ———————— 需要正规操作，如果以后有时间再修改
% I路编码信息
data_code = [data_i zeros(1,padding)];                          % apply CRC
channel_data_i = convenc(data_code, trellis);                   % apply convlution code
ch_i_map = channel_data_i*2-1;
% Q路编码信息
channel_data_q=[1 1 1 1 1 0 1 1 1 1,1 0 0 1 1 0 1 1 1 1,1 0 1 1 0 1 1 1 1 1,...
       1 0 0 1 0 0 1 1 1 1,1 1 0 1 0 1 1 1 1 1,1 1 1 1 1 0 1 1 1 1,...
       1 1 1 1 0 0 1 1 1 1,1 1 0 1 0 0 1 1 1 1,1 0 1 1 1 0 1 1 1 1,...
       1 1 1 1 1 1 1 1 1 1,1 0 1 1 0 1 1 1 1 1,1 1 0 1 1 1 1 1 1 1,...
       1 1 0 1 0 0 1 1 1 1,1 0 0 1 1 1 1 1 1 1,1 0 0 1 1 1 1 1 1 1];
ch_q_map = channel_data_q*2-1;

if length(channel_data_i)~=sim_channel_i || length(channel_data_q)~=sim_channel_q
    disp("error:信道编码长度不匹配");
    return;
end

%% 加扩
% I路加扩
s_spread_i = [];
for i = 1 : length(ch_i_map)
    s_spread_i = [s_spread_i ch_i_map(i) * ovsf_codes_i];
end

% Q路加扩
s_spread_q = [];
for i = 1 : length(ch_q_map)
    s_spread_q = [s_spread_q ch_q_map(i)*ovsf_codes_q];
end

if length(s_spread_i)~= sim_channel_i*spread_factor_i || length(s_spread_q)~= sim_channel_q*spread_factor_q
    disp("error:扩频操作长度不匹配");
    return;
end

% 加权
bd=1.0;bc=1.0;
spread_s = bd*s_spread_i + j*bc*s_spread_q;  

%% 加扰
load ./datafile/gold_code.mat;                                  % 长度为38400的gold码
gold = gold_code;
scramble_s = spread_s.*gold; 


%% 成型滤波
s_real = real(scramble_s);
s_imge = imag(scramble_s);

rolloff_factor = 0.5;                                           % 滚降因子
rcos_fir = rcosdesign(rolloff_factor,2*sps,sps);                % 默认是根升余弦滤波器,'sqrt'

% 插值
for k=1:length(s_real)
    up_sps_s_real(1+sps*(k-1)) = s_real(k);
    up_sps_s_real(2+sps*(k-1):sps*k) = zeros(1,sps-1);
end

for k=1:length(s_imge)
    up_sps_s_imag(1+sps*(k-1)) = s_imge(k);
    up_sps_s_imag(2+sps*(k-1):sps*k) = zeros(1,sps-1);
end

% 滚降滤波
rcos_s_real = filter(rcos_fir,1,up_sps_s_real);
rcos_s_imag = filter(rcos_fir,1,up_sps_s_imag);

rcos_scramble_s = rcos_s_real + j*rcos_s_imag;


%% 加噪声段
pre_noise = zeros(1,384000);
signal_burst = [pre_noise rcos_scramble_s];
%% 信道
signal_burst_noise = awgn(signal_burst,snr,'measured');
send_signal = signal_burst_noise./max(abs(signal_burst_noise));

end
