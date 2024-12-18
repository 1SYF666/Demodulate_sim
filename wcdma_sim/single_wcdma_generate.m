%{
    ���������
            fream_t--֡ʱ����fs--�����ʣ�Rb--��Ϣ���ʣ�Rb_I--I·�ŵ����ʣ�
            spread_factor_i--I·��Ƶ�볤��sf_num--��Ƶ���ţ�SNR--����ȣ�
            Q·��ϢΪ�����ŵ���Q·��Ϣ����������15K���Լ���Ƶ�볤Ϊ256���ǹ̶��ġ�
    ���������
            send_signal--���������źţ�gold--���룻ovsf_codes_i--I·��Ƶ�룻
            ovsf_codes_q--Q·��Ƶ�룻channel_data_i--�ŵ������i·���ݣ�
            channel_data_q--�ŵ������q·���ݣ�data_code -- CRC����Ϣ���ݣ�
            data_i -- ��Ϣ����
    I·��Ϣ���ʣ�600  2400 32000 64000 128000 384000 
%}
function [send_signal,gold,ovsf_codes_i,ovsf_codes_q,...
    channel_data_i,channel_data_q,data_code,data_i]...
    = single_wcdma_generate(fream_t,Rb,Rb_I,spread_factor_i,sf_num,snr,sps)  
%% ��������  ������������  ��Ҫ���ŵ����룬�������ʱ��Ļ�
j = sqrt(-1);
% I·��������  
if Rb==600
    padding = 28+16;                
    trellis = poly2trellis(9, [771 663 557]);                   %��3��1��9���������
elseif Rb==2400
    padding = 10+16;                
    trellis = poly2trellis(9, [771 663 557]);                   %��3��1��9���������
elseif Rb==32000
    padding=64+16;
    trellis = poly2trellis(9, [771 663 557]);                   %��3��1��9���������
elseif Rb==64000
    padding=144+16;
    trellis = poly2trellis(9, [771 663 557]);                   %��3��1��9���������
elseif Rb==128000
   padding=304+16;
   trellis = poly2trellis(9, [771 663 557]);                    %��3��1��9���������
elseif Rb==384000
    if spread_factor_i==4
        trellis = poly2trellis(9, [561 753]);                     %��2��1��9���������
        padding=944+16;
    elseif spread_factor_i==2
        trellis = poly2trellis(9, [561 753]);
        padding=5760;
    end
end
sim_base_i = fream_t*Rb;                                        % ��Ϣ���ظ���
sim_channel_i = fream_t*Rb_I;                                   % �ŵ����ظ���
ovsf_codes_i = ovsf_generator(spread_factor_i,sf_num);          % I·��Ƶ��

% Q·��������
Rb_Q = 15000;                                                   % Q·�ŵ�����
spread_factor_q = 256;              
sim_channel_q = fream_t*Rb_Q;
ovsf_codes_q = ovsf_generator(spread_factor_q,1);               % Q·��Ƶ��

%% ��Ϣ����
data_i = randi([0 1],1,sim_base_i);                             % in 10 ms, 15 slots

%% �ŵ����� ���������������� ��Ҫ�������������Ժ���ʱ�����޸�
% I·������Ϣ
data_code = [data_i zeros(1,padding)];                          % apply CRC
channel_data_i = convenc(data_code, trellis);                   % apply convlution code
ch_i_map = channel_data_i*2-1;
% Q·������Ϣ
channel_data_q=[1 1 1 1 1 0 1 1 1 1,1 0 0 1 1 0 1 1 1 1,1 0 1 1 0 1 1 1 1 1,...
       1 0 0 1 0 0 1 1 1 1,1 1 0 1 0 1 1 1 1 1,1 1 1 1 1 0 1 1 1 1,...
       1 1 1 1 0 0 1 1 1 1,1 1 0 1 0 0 1 1 1 1,1 0 1 1 1 0 1 1 1 1,...
       1 1 1 1 1 1 1 1 1 1,1 0 1 1 0 1 1 1 1 1,1 1 0 1 1 1 1 1 1 1,...
       1 1 0 1 0 0 1 1 1 1,1 0 0 1 1 1 1 1 1 1,1 0 0 1 1 1 1 1 1 1];
ch_q_map = channel_data_q*2-1;

if length(channel_data_i)~=sim_channel_i || length(channel_data_q)~=sim_channel_q
    disp("error:�ŵ����볤�Ȳ�ƥ��");
    return;
end

%% ����
% I·����
s_spread_i = [];
for i = 1 : length(ch_i_map)
    s_spread_i = [s_spread_i ch_i_map(i) * ovsf_codes_i];
end

% Q·����
s_spread_q = [];
for i = 1 : length(ch_q_map)
    s_spread_q = [s_spread_q ch_q_map(i)*ovsf_codes_q];
end

if length(s_spread_i)~= sim_channel_i*spread_factor_i || length(s_spread_q)~= sim_channel_q*spread_factor_q
    disp("error:��Ƶ�������Ȳ�ƥ��");
    return;
end

% ��Ȩ
bd=1.0;bc=1.0;
spread_s = bd*s_spread_i + j*bc*s_spread_q;  

%% ����
load ./datafile/gold_code.mat;                                  % ����Ϊ38400��gold��
gold = gold_code;
scramble_s = spread_s.*gold; 


%% �����˲�
s_real = real(scramble_s);
s_imge = imag(scramble_s);

rolloff_factor = 0.5;                                           % ��������
rcos_fir = rcosdesign(rolloff_factor,2*sps,sps);                % Ĭ���Ǹ��������˲���,'sqrt'

% ��ֵ
for k=1:length(s_real)
    up_sps_s_real(1+sps*(k-1)) = s_real(k);
    up_sps_s_real(2+sps*(k-1):sps*k) = zeros(1,sps-1);
end

for k=1:length(s_imge)
    up_sps_s_imag(1+sps*(k-1)) = s_imge(k);
    up_sps_s_imag(2+sps*(k-1):sps*k) = zeros(1,sps-1);
end

% �����˲�
rcos_s_real = filter(rcos_fir,1,up_sps_s_real);
rcos_s_imag = filter(rcos_fir,1,up_sps_s_imag);

rcos_scramble_s = rcos_s_real + j*rcos_s_imag;


%% ��������
pre_noise = zeros(1,384000);
signal_burst = [pre_noise rcos_scramble_s];
%% �ŵ�
signal_burst_noise = awgn(signal_burst,snr,'measured');
send_signal = signal_burst_noise./max(abs(signal_burst_noise));

end
