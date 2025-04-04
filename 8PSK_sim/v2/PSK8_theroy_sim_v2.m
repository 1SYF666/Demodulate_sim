%% ************* 8PSK simulation ************* %%
%% ***** data:20240929 authoor:ShenYifu ****  %%
%% 参考资料《EPDT空中接口物理层和数据链路层》中物理层8PSK调制
%% 参考资料《GMSK、8PSK调制算法研究及ASIC实现》_肖妍妍
%% 参考资料《GSM/EDGE 与 8PSK调制信号性能比较》_罗宏杰
%% 参考资料《GSM/EDGE基带关键模块的设计与实现》_常培峰
%% note:解旋之后直接判决
clc;clear;
close all;
j = sqrt(-1);
%% 基本参数
% Rb = 50e3;
Rb = 1625*1e3/6;
sps = 4;
fs = sps * Rb;
Ts = 1 / fs;
Tb = 1 / Rb;
EbN0 = 16;
DEBUG = 0;
%% 8PSK调制
base_bit = load('C:\Users\PC\Desktop\postgraduate\others\东方通信\8psk_dongfang\psk8_signal\50K.txt');
[s,si,base_symbol] = pskmod(base_bit);
%% 成型滤波
s_real = real(s);
s_imge = imag(s);

C_0 = pulsefilter(Tb,Ts*1.5,5);

for i = 1 : length(s_real)*sps
    y_temp = 0;
    for  k = 1 : length(s_real)
        index = i-k*(sps)+2.5*sps;
        if index < 1||index>length(C_0)
            continue;
        end
        c0temp = C_0(index);
        y_temp = y_temp + s(k)*c0temp;
    end
    s_rcos(i) = y_temp;
end

%% 同步数据
synciq1 =  si(50+4+208+1 : 50+4+208+24);
synciq2 =  si(50+4+208+24+208+1 : 50+4+208+24+208+24);
synciq = [synciq1 zeros(1,208) synciq2];
Gp_start = s(1:50);
Tb_start = s(50+1:50+4);
sync1 =  s(50+4+208+1 : 50+4+208+24);
sync2 =  s(50+4+208+24+208+1 : 50+4+208+24+208+24);
Tb_end = s(50+4+208+24+208+24+208+1:50+4+208+24+208+24+208+4);
Gp_end = s(50+4+208+24+208+24+208+4+1:50+4+208+24+208+24+208+4+20);
sync = [sync1 zeros(1,208) sync2];

%% 信道
noise = 1;
send = repmat(s_rcos,1,10);   
send_signal = awgnself(send,noise,EbN0,sps);
%% 8PSK解调
receive = send_signal;

%% 匹配滤波
receive_rcos = receive;

%% 粗同步
[headframe,tailframe] = coarse_synchronisation(receive_rcos,sync,sps);

%% 分帧处理
posA = [];
maxPosA = [];
posArraryA = [];
framestart = 208+4+50; % 根据8PSK帧结构而定
for i = 1 : 5
    framesignal = [];
    framesignal = receive_rcos(headframe(i):tailframe(i));

    %% 符号同步
    maxC = 0; pos=0; maxR=0;
    meanA = mean(abs(framesignal));
    maxC = zeros(1,sps);
    posArray = zeros(1,sps);
    for ddcIdx=1:sps   % 1:18
        a = xcorr(framesignal(ddcIdx:sps:end),sync);  % 仅需滑动1~3个样点，主要是计算相位值maxC
        % figure;plot(abs(a));title("同步自相关示意图",ddcIdx);
        [maxV,p]=max(abs(a));
        maxR = maxV;
        pos = (p-length(framesignal(ddcIdx:sps:end))-framestart)*sps+ddcIdx;
        maxC(ddcIdx) = a(p);
        posArray(ddcIdx) = pos;
    end
    % 找最大峰值ddc位置
    [~,maxPos] = max(abs(maxC));
    maxPosA = [maxPosA,maxPos]; % 修改,记录maxPosA的值
    posArraryA = [posArraryA,posArray];
    pos = posArray(maxPos);
    posA = [posA,pos];  % 修改,记录posA的值
   
    if pos<1
        continue;  
    end
    s_synchronization = framesignal(pos:sps:end)*maxC(maxPos)';

    %% 解旋
    rotsignal = s_synchronization;
    if framestart == 208+4+50
        k = 1 : 750;
        ppp = exp(-1i*3*pi/8*k);
        for k = 1 : length(ppp)
            s_derot(k) = rotsignal(k).*ppp(k);
        end
    end


    %% 解映射
    y = desymbolmap(s_derot);
    c_symbol_N = y;
    [deoutbit] = debitmap(y);
    %% 误码率
    errorsymbol = sum(c_symbol_N~= base_symbol); 
    fprintf("第 %d 帧：误符号数：%d ,误符号率为：%.6f ; ",i,errorsymbol,errorsymbol/length(c_symbol_N));
    % figure;stem(c_symbol_N(4:end-3));hold on;stem(base_symbol(4:end-3));title("解调输出符号对比图");
    % figure;plot(y);hold on ;plot(0.95*base_symbol);title("解调输出符号对比示意图");
    errorbit = sum(deoutbit~=base_bit); 
    fprintf("误比特数：%d ,误比特率为：%.6f\n",errorbit,errorbit/length(base_bit));
    % figure;stem(deoutbit(1:end));hold on;stem(base_bit(3*3+1:end-3*3));title("解调输出比特对比图");
    % figure;plot(deoutbit(1:end));hold on;plot(0.98*base_bit(3*3+1:end-3*3));title("解调输出比特对比图");
end
figuredebug;
