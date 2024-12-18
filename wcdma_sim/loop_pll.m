%% 锁相环去频偏
% 对解扰信号进行积分操作，选择256的原因：
% 首先：I路扩频码为长度(256、128、64、32等)，每一个码内+1和-1个数都是一样的，
% 以至于对一个符号或者多个符号进行累加的话相当于 这一个符号或多个符号都变成0了
% 其次：Q路扩频码长为256，且全为1，这个长度累加起来相当于还原这个符号了

% 第一次积分操作--目的：解Q路数据
sum_length=spread_factor_I*sps;
st=1;
sum_number=floor(length(Vo_des)/(sum_length));
for i=st:sum_number
    Vo_sum(i)=sum(Vo_des((i-1)*sum_length+1:i*sum_length))/sum_length; %除以sum_length相当于归一化吧？暂时想不出合理的解释
    Vo_dep_A_sum(i)=sum(Vo_dep_A((i-1)*sum_length+1:i*sum_length))/sum_length; %目的是为了解I路数据

end
figure;scatter(real(Vo_sum),imag(Vo_sum)); title("进锁相环前星座图");
figure;subplot(2,1,1);plot(real(Vo_sum));title("进锁相环前I路效果图");subplot(2,1,2);plot(imag(Vo_sum));title("进锁相环前Q路效果图");
    
Simulation_Length=length(Vo_sum);
Signal_Channel=zeros(1,Simulation_Length);
Signal_PLL=zeros(1,Simulation_Length);NCO_Phase =zeros(1,Simulation_Length);
Discriminator_Out=zeros(1,Simulation_Length);Freq_Control=zeros(1,Simulation_Length);
PLL_Phase_Part=zeros(1,Simulation_Length);PLL_Freq_Part=zeros(1,Simulation_Length);

Ko=1;Kd=1;K=Ko*Kd;
sigma=0.707;
symbol_rate=sps*3.84e6/sum_length;  % 采样率除以积分点数表示进锁相环时的符号速率
BL=0.15*symbol_rate;                  %频率0.0254--
Wn=8*sigma*BL/(1+4*sigma^2);     T_nco=1/symbol_rate;
C11=(2*sigma*Wn*T_nco)/(K);      C22=((T_nco*Wn)^2)/(K);

Vi_PLL=Vo_sum;                       %解扰信号输入锁相环
Vi_PLL_real=real(Vi_PLL);  Vi_PLL_imag=imag(Vi_PLL);

for i=2:Simulation_Length
    Signal_Channel(i)=Vi_PLL_real(i)+1i*Vi_PLL_imag(i);
    Signal_PLL_D(i)=Signal_Channel(i)*exp(-1i*(mod(NCO_Phase(i-1),2*pi)));    
    Vo_PLL_real(i)=real(Signal_PLL_D(i));  %科斯塔斯环输出
    Vo_PLL_imag(i)=imag(Signal_PLL_D(i));
    Discriminator_Out(i)=sign(Vo_PLL_real(i))*Vo_PLL_imag(i)-sign(Vo_PLL_imag(i))*Vo_PLL_real(i);
    PLL_Phase_Part(i)=Discriminator_Out(i)*C11;
    PLL_Freq_Part(i)=Discriminator_Out(i)*C22+PLL_Freq_Part(i-1);
    Freq_Control(i)=PLL_Phase_Part(i)+PLL_Freq_Part(i);
    NCO_Phase(i)=NCO_Phase(i-1)+Freq_Control(i)*2*pi;
    
    Vo_dep(i)=Vo_dep_A_sum(i).*exp(-1i*(mod(NCO_Phase(i-1),2*pi)));  %保存相位变化值，进行下面对I路数据的解扩
end

figure;plot(PLL_Freq_Part);title("载波同步频率跟踪曲线");

% 由于Q路扩频码长度为256，而第一次积分以spread_factor_I长度积分,
% 那么锁相环一次输出并代表q路一个符号，而是代表的是spread_factor_I/256个符号
% 故再需要对锁相环输出一次累加长度为256/spread_factor_I，此时的一次输出才能表示q路的一个符号
% 第二次积分操作--对第一次的补充
sum_length2=256/spread_factor_I;
st2=1;
sum_number2=floor(length(Vo_PLL_real)/(sum_length2));   % Vo_dep_A---解扰后信号又乘以I路的扩频码，为了使Q路累加为0
for i=st2:sum_number2
    Vo_sum2(i)=sum(Vo_PLL_real((i-1)*sum_length2+1:i*sum_length2))/sum_length2; %除以sum_length2相当于求出一个符号的值了
end
figure;plot(Vo_sum2);title("Q路数据输出结果");
figure;plot(2*Vo_sum2);hold on;plot(1.2*despread_q);title("解扩后与加扩前Q路数据对比");

% 画出I路数据
figure;plot(real(Vo_dep));title("I路数据输出结果");
figure;plot(2*real(Vo_dep));hold on;plot(1.2*despread_i);title("解扩后与加扩前I路数据对比");
