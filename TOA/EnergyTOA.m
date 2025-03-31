clc;clear;
close all;
j = sqrt(-1);
%% ��������
K = 1e3;             % ��λ KHz
M = 1e6;             % ��λ MHz
Rb = 1*K;            % ������
fs = 4 * Rb;      % ������
time = 5;            % ����ʱ��
symbolnum = Rb*time; % ��Ԫ����
fc = 1500*M;       
SNR_dB = 15;        % �����(��λ dB);
start=30000;
finish=30000;
T = 1/fs;           % ������� (s)
starttime=start*T;
t = 0:T:1-T;        % ʱ������ (1s)
% T_b = 0.1;          % �����ۻ����� (s)
% W = round(T_b / T); % ÿ��������Ĳ�������
W = 32;      %�������������   
W_slip = 32;   %�ƶ�����
threadhold1 = 28;   %�ƶ�������ż�
toa_est = [];  %��ʼ��
toe_est = [];  %��ֹ��
startpoint = [];
finishpoint = [];
Esum = [];
Ysum = [];
start_counter = 0;  % ��ʼ�������
finish_counter = 0; % ��ֹ�������

%% BPSK�ź�����
% PCM����
pcm_bits = randi([0 1],1,symbolnum);  %���ɳ���Ϊsymbolnum�������randi������������(pcm_bits)��ֵΪ0��1
pcm_symbols = 2 * pcm_bits - 1;   % ��0 --> -1; 1 --> 1   %pcm_symbols,ӳ����BPSK����
% BPSK ����
sps = fs / Rb;                    % spsÿ���Ų�����
rcos_pcm_s1 = repmat(pcm_symbols,sps,1);     %repmat��ÿ�������ظ�sps��
rcos_pcm_s1(2:sps,:) = 0;      %���ڶ��е���sps�е�ֵ��Ϊ0
rcos_pcm_s = rcos_pcm_s1(:)';   %������չƽΪ������
rolloff_factor = 0.7;       % �������ӣ�ȡֵ��Χ[0,1]�������˲����Ĵ�����԰�˥��
rcos_fir = rcosdesign(rolloff_factor,2*sps,sps); % Ĭ���Ǹ��������˲���,'sqrt'
%2*sps���˲����Ľ���������Ϊ2*sps+1��rcos_fir�����������˲���ϵ��
rcos_ds_real = filter(rcos_fir,1,rcos_pcm_s);  %1���˲����ķ�ĸϵ����rcos_pcm_real��������κ���ź�

%% ��Դ�ź�ģ��
S1 = rcos_ds_real+j*rcos_ds_real;  %S1Ϊ���źţ�ʵ��+�鲿
%fa = exp(-j*2*pi*fc*d*sin(theta)/c);
%disp(['S1 �ĳ��ȣ�', num2str(length(S1))]);  % ��ӡ S1 �ĳ���
Y = [zeros(1,start) S1  zeros(1,finish) S1 zeros(1,finish) S1 zeros(1,finish)];   %���ɷ����źţ�ǰ���3W���㣨����䣩
%figure;plot(real(Y));   %����������Y��ʵ��

%% �ŵ�
SNR1_dB = SNR_dB-log(sps)/log(10)*10;  
%SNR_dBԭʼ����ȣ�log(sps)/log(10)*10��spsת��Ϊ�ֱ���SNR1_dB������������
%��ʽ��SNR1_dB=SNR_dB?10?log10��sps��
P_sig = mean( abs(S1(:)).^2 );   
%S1(:)��S1չƽΪ��������abs(S1(:)).^2�����ź�ÿ��������Ĺ��ʣ�mean����ƽ������
P_noise = P_sig / (10^(SNR1_dB/10));
%10^(SNR1_dB/10)������ȴӷֱ�ת��Ϊ����ֵ�����幫ʽ����������ȼ�����������
if 0    
    % �����ʵֵ����:
    noise_matrix = sqrt(P_noise) * randn(1, length(Y));
    %randn(1,length(Y)):���ɷ��ӱ�׼��̬�ֲ��������������Ϊ1��
    ...........................................................++++++++++
else
    % ����Ǹ�ֵ����:
    noise_matrix = sqrt(P_noise/2) * (randn(1, length(Y)) + 1i*randn(1, length(Y)));
    %sqrt(P_noise/2):�����������ȣ�����ΪP_noise;1i:������λ
end

Y_noise = Y + noise_matrix;
figure;plot(real(Y));
title('�����źţ���������');

%% ��8192�ֶδ����ź�
segment_length=8192;
sign_length=length(Y_noise);
num_segments = floor(sign_length / segment_length)+1;

%% 
startflag=0;
for k = 1:num_segments
    % ��ȡ��ǰ��
    start_index = (k - 1) * segment_length + 1;  %��ǰ���źŵ���ʼ��
    end_index = min(k * segment_length, sign_length);;    %��ǰ���źŵ���ֹ��
    segment = Y_noise(start_index:end_index);   %��ȡָ��������ź�
    segment_Y = Y(start_index:end_index);   %��ȡ�������������ź�
    %% ������ⷨ
    % ƽ���ʼ첨
    segment_squared = abs(segment).^2;

    % �����ۻ�
    num_blocks = floor(length(segment_squared) / W);  % ����������
    E = zeros(1, num_blocks);  % ������
    for lambda = 1:num_blocks
        start_idy = (lambda-1)*W + 1;
        end_idy = lambda*W;
        E(lambda) = sum(segment_squared(start_idy:end_idy));  % �����ۻ�
    end
    % ��չ�����ۻ��ź�
    ones1 = ones(W,1);
    E_spread = ones1*E;
    E_spread2 = E_spread(:).';

    % ���� E_spread2 �ĳ����� Y ��ͬ
    if length(E_spread2) < length(segment)
        E_spread2 = [E_spread2, zeros(1, length(segment) - length(E_spread2))];  % ����
    else
        E_spread2 = E_spread2(1:length(segment));  % �ض�
    end
    Esum = [Esum,  E_spread2];

    %% �����ź�Y��E

    % ������ԭʼ�ź�
    downsample_factor = 1;  % ����������
    Y_downsampled = segment_Y(1:downsample_factor:end);  % ����� downsample_factor ��ȡ������

    % �����������ۻ��ź�
    E_spread2_downsampled = E_spread2(1:downsample_factor:end);



    %% ���ù�һ������
    Z_max = max(E);  % �������е����ֵ
    Z_min = min(E);  % �������е���Сֵ
    gamma_norm = 0.5;  % ��һ�����ޱ���
    gamma_facc = gamma_norm * (Z_max - Z_min) + Z_min;  % �̶�����
    %����tʱ����Y�ĳ��ȶ�̬����
    t = (0:length(segment)-1) * T;  % ���� Y �ĳ�������ʱ������

    %% �ƶ������
    Wtemp = zeros(1,W_slip);
    for i =  1 : num_blocks-W_slip
        Wtemp =  E(i:i+W_slip-1);  %��������ֵ���ƶ���
        Wtemp = Wtemp>gamma_facc;
        gamma_slip = sum(Wtemp);  %�����������
        %Ѱ����ʼ��
        if  gamma_slip >  threadhold1 && startflag == 0%����������ֵ�Ƚ�
            start_counter = start_counter + 1;
            toa_est = i ;  %��λ��ʼ��  
            startflag = 1;  %֤���Ѿ��ҵ���ʼ�㣬��ʼѰ�ҽ�����
            disp(['�ҵ���',num2str(start_counter),'����ʼ��']);
            startpoint = [startpoint, i*W + (k-1)*segment_length];  %��ʼ������
        end
        %Ѱ����ֹ��
        if gamma_slip < W_slip-threadhold1 && startflag == 1
           finish_counter = finish_counter + 1;  % ��ֹ���������1
           toe_est = i;  %��λ��ֹ��
           startflag = 0;   %֤���Ѿ��ҵ���ֹ�㣬��ʼѰ����ʼ��
           disp(['�ҵ��� ', num2str(finish_counter), ' ����ֹ��']);
           finishpoint = [finishpoint, i*W + (k-1)*segment_length];  %��ֹ������
        end
    end

end
    %% ������������
    figure;
    plot(Esum);  % ��һ�������ۻ��ź�
    hold on;
    plot(real(Y_noise));  % ��һ��ԭʼ�ź�
    title('�����ۻ��ź���ԭʼ�ź� ');
    xlabel('������');
    ylabel('��һ������');
    legend('�����ۻ��ź�', 'ԭʼ�ź�');
%% ���
    figure;plot(real(Y_noise));  
    hold on;
    plot(startpoint, real(Y_noise(toa_est)), 'ro');
    plot(finishpoint,real(Y_noise(toa_est)), 'gs');
    title('�ź���ʼ��ֹ����');
    xlabel('������');
    legend('�����ź�', '���Ƶ� TOA', '���Ƶ� TOE');