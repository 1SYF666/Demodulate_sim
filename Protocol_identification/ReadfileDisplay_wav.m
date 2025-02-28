%% ********* Read file and frequency_time figure ********* %%
%% ********* data:20250228 authoor:ShenYifu ****  %%
%% ***** I类信号文件  ***** %
%{
    I类信号文件中含有格式：.spl ; .wav ; std 格式
    其中 .spl格式： 是复信号，int16
    其中 .std格式： 是复信号，int32
    其中 .wav格式： 实信号或者复信号
    wav文件名格式如下：
    IF_2019-06-11_16-36-42_294175000Hz_120_kHz.wav
%}

%%
clc;clear;
close all;
%% 参数设置
K = 1e3;           % 单位1KHz
%% 读取.wav格式文件
folderPath1 = 'D:\MODIFI_CODE\协议信号研究\I'; 
files1 = dir(fullfile(folderPath1, '*.wav')); 

folderPath = 'D:\MODIFI_CODE\协议信号研究\I类信号文件截取';
files = dir(fullfile(folderPath, '*.wav')); 

%% 保存图片地址
folderPath2 = 'D:\MODIFI_CODE\协议信号研究\I1\wave';  % 修改为你的文件夹路径
if ~exist(folderPath2, 'dir')
    mkdir(folderPath2);
end
folderPath3 = 'D:\MODIFI_CODE\协议信号研究\I1\timefrequency';  
if ~exist(folderPath3, 'dir')
    mkdir(folderPath3);
end

%%
for j = 2:length(files1)
    % 获取文件名
    fileName1 = files1(j).name;
    targetPrefix = fileName1; 
    [~, targetPrefix1, ~] = fileparts(targetPrefix);

    % 获取信号带宽
    pattern = '(?<=_)(\d+)(?=_kHz)';% 使用正则表达式提取kHz和_之间的数字部分
    number = regexp(targetPrefix1, pattern, 'match');% 输出提取到的数字
    if ~isempty(number)
        band_extract = str2double(number{1}) * K;
        disp(['提取的带宽是: ', num2str(band_extract)]);
    else
        disp('未找到符合条件的数字');
        break;
    end
    % 读取wav格式信息
    filePath1 = fullfile(files1(j).folder, files1(j).name);
    fid1 = fopen(filePath1, 'rb');
    fseek(fid1,22,'bof');
    signal_num_channels =fread(fid1,1,'uint16'); % 单(双)通道
    fseek(fid1,34,'bof');
    signal_bit =fread(fid1,1,'uint16');          % 数据格式
    fseek(fid1,24,'bof');
    signal_fs=fread(fid1,1,'int32');             % 采样率
    fclose(fid1);

    % 删除文件夹中现有的所有图像文件（根据扩展名过滤）
    delete(fullfile(folderPath2, '*.png'));  % 删除文件夹中的所有.png文件
    delete(fullfile(folderPath2, '*.jpg'));  % 如果有.jpg文件，也可以删除
    delete(fullfile(folderPath2, '*.fig'));  % 如果有.fig文件，也可以删除
    delete(fullfile(folderPath3, '*.png'));
    delete(fullfile(folderPath3, '*.jpg'));
    delete(fullfile(folderPath3, '*.fig'));
    
    figurenum = 0;
    % 遍历文件夹中的每个文件
    for i = 1:length(files)
        fileName = files(i).name;
        fileContent = fileName;

        if strncmp(fileContent, targetPrefix1, length(targetPrefix1))
            disp('匹配成功，执行后续操作');
            disp(['Processing file: ', fileName]);
            filePath = fullfile(folderPath, fileName);
            fid = fopen(filePath, 'rb');
            
            % 读取wav数据部分
            fseek(fid,44,'bof'); 
            if signal_bit==8
                data=fread(fid,inf,'int8');
            elseif signal_bit==16
                data=fread(fid,inf,'int16');
            elseif signal_bit==32
                data=fread(fid,inf,'int32');
            end
            fclose(fid);
            
            % 单(双)通道
            if signal_num_channels==1
                signal=data.';
            elseif signal_num_channels==2
                signal=data(1:2:end-1)+1i*data(2:2:end);
                signal=signal.';
            end
            
            % 时频图参数设置
            N = length(signal);
            window = floor(N / 10);  % 设置窗口长度为信号长度的1/10
            noverlap = floor(window / 2);  % 设置重叠样本数为窗口的一半
            nfft = 2^nextpow2(window);  % 设置FFT点数为2的下一个最接近的整数次幂
            if window > N
                window = N;
                noverlap = N - 1;
            end

            % 绘图
            figurenum = figurenum + 1;
            if signal_num_channels == 1
                figure('Visible','off');
                plot(signal);title("波形图");
            elseif signal_num_channels == 2
                figure('Visible','off');
                subplot(2,1,1);plot(real(signal));title("实部和虚部波形图");subplot(2,1,2);plot(imag(signal));
                signal = signal.*exp(-1i*2*pi*(signal_fs/2)*(1:length(signal))/signal_fs);
            end
            filesave = fullfile(folderPath2, ['figure' num2str(figurenum) '.png']);
            saveas(gcf, filesave);  % 保存图形为png格式

            % [s, f, t] = spectrogram(sdata,window,noverlap,nfft,fs_extract);
            figure('Visible','off'); 
            spectrogram(signal,window,noverlap,nfft,signal_fs,'yaxis');title('Time-Frequency Plot');
            filesave2 = fullfile(folderPath3, ['figure' num2str(figurenum) '.png']);
            title("");
            axis off;
            set(gca, 'XTick', [], 'YTick', []);% 去掉坐标轴刻度
            colorbar off; % 去掉颜色条
            set(gca, 'Position', [0, 0, 1, 1]);  % 调整图形的边距，去除白边,这行会移除额外的空白区域
            
            saveas(gcf, filesave2);  % 保存图形为png格式

            close all;
        end
    end
end