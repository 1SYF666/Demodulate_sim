Simulation_Length=length(Vo_sum);
Signal_Channel=zeros(1,Simulation_Length);
Signal_PLL=zeros(1,Simulation_Length);NCO_Phase =zeros(1,Simulation_Length);
Discriminator_Out=zeros(1,Simulation_Length);Freq_Control=zeros(1,Simulation_Length);
PLL_Phase_Part=zeros(1,Simulation_Length);PLL_Freq_Part=zeros(1,Simulation_Length);

Ko=1;Kd=1;K=Ko*Kd;
sigma=0.707;
symbol_rate=looprb  ;% 采样率除以积分点数表示进锁相环时的符号速率
BL=loopcoef*symbol_rate;                  %频率0.0254--
Wn=8*sigma*BL/(1+4*sigma^2);     T_nco=1/symbol_rate;
C11=(2*sigma*Wn*T_nco)/(K);      C22=((T_nco*Wn)^2)/(K);

