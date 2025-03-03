%% ********* Read file and frequency_time figure ********* %%
%% ********* data:20250228 authoor:ShenYifu ****  %%
%% ***** I类信号文件  ***** %
%{
    读取std转换后文件
    std文件名格式如下：
    T610X_DDCRX1_2019-12-04-05-29-30_F157.625MHz_BW640kHz_Fs800Ksps.std
%}

%%

clc;clear;
close all;

%% 参数设置
K = 1e3;           % 单位1KHz
%% 读取.spl格式文件
folderPath1 = 'D:\MODIFI_CODE\协议信号研究\B'; % % 指定文件夹路径
files1 = dir(fullfile(folderPath1, '*.std')); % 获取该文件夹下所有文件的列表
folderPath = 'D:\MODIFI_CODE\协议信号研究\B类信号文件截取'; % % 指定文件夹路径
files = dir(fullfile(folderPath, '*.dat')); % 获取该文件夹下所有文件的列表

%% 保存图片地址
folderPath2 = 'D:\MODIFI_CODE\协议信号研究\B1\wave';  % 修改为你的文件夹路径
if ~exist(folderPath2, 'dir')
    mkdir(folderPath2);
end
folderPath3 = 'D:\MODIFI_CODE\协议信号研究\B1\timefrequency';  
if ~exist(folderPath3, 'dir')
    mkdir(folderPath3);
end


%%
for j = 1:length(files1)
% 获取文件名
    fileName1 = files1(j).name;
    targetPrefix = fileName1;  % 要匹配的前缀字符串
    [~, targetPrefix1, ~] = fileparts(targetPrefix);% 使用fileparts提取文件名和扩展名

    % 获取采样率
    % pattern = '(?<=_)(\d+)(?=sps)';% 使用正则表达式提取sps和_之间的数字部分
    % number = regexp(targetPrefix1, pattern, 'match');% 输出提取到的数字
    % if ~isempty(number)
    %     signal_fs = str2double(number{1});
    %     disp(['提取的采样率是: ', num2str(signal_fs)]);
    % else
    %     disp('未找到符合条件的数字');
    %     break;
    % end
    signal_fs = 80e3;
    % close all; % 清除图片
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

        % 检查文件内容是否以指定的前缀开头
        if strncmp(fileContent, targetPrefix1, length(targetPrefix1))
            disp('匹配成功，执行后续操作');

            filePath = fullfile(folderPath, fileName);

            % 读取文件内容（根据文件格式调整读取方法）
            fileID = fopen(filePath, 'r');
            fileData = fread(fileID, 'int16')';
            fclose(fileID);

            % 对文件数据进行处理
            disp(['Processing file: ', fileName]);

            idata = fileData(1:2:end-1);            % 复信号
            qdata = fileData(2:2:end);
            signal = idata+1i*qdata;
            N = length(signal);
            window = floor(N / 10);  % 设置窗口长度为信号长度的1/10
            noverlap = floor(window / 2);  % 设置重叠样本数为窗口的一半
            nfft = 2^nextpow2(window);  % 设置FFT点数为2的下一个最接近的整数次幂

            if window > N
                window = N;            % 确保窗口长度不大于信号长度
                noverlap = N - 1;
            end

            % 绘图
            figurenum = figurenum + 1;

            figure('Visible','off');
            subplot(2,1,1);plot(real(signal));title("实部和虚部波形图");subplot(2,1,2);plot(imag(signal));
            signal = signal.*exp(-1i*2*pi*(signal_fs/2)*(1:length(signal))/signal_fs);

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