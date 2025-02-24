%% ********* Read file and frequency_time figure ********* %%
%% ********* data:20250224 authoor:ShenYifu ****  %%
%% ***** I类信号文件  ***** %
%{
    I类信号文件中含有格式：.spl ; .wav ; std 格式
    其中 .spl格式： 是复信号，int16
    其中 .std格式： 是复信号，int32
    其中 .wav格式： (暂时未知)
%}

%%
clc;clear;
close all;
%% 读取.spl格式文件
folderPath1 = 'D:\MODIFI_CODE\协议信号研究\I'; % % 指定文件夹路径
files1 = dir(fullfile(folderPath1, '*.spl')); % 获取该文件夹下所有文件的列表

folderPath = 'D:\MODIFI_CODE\协议信号研究\I类信号文件截取'; % % 指定文件夹路径
files = dir(fullfile(folderPath, '*.pcm')); % 获取该文件夹下所有文件的列表
for j = 1:length(files1)
    % 获取文件名
    fileName1 = files1(j).name;
    targetPrefix = fileName1;  % 要匹配的前缀字符串
    [~, targetPrefix1, ~] = fileparts(targetPrefix);% 使用fileparts提取文件名和扩展名

    % 获取采样率
    pattern = '(?<=_)(\d+)(?=sps)';% 使用正则表达式提取sps和_之间的数字部分
    number = regexp(targetPrefix1, pattern, 'match');% 输出提取到的数字
    if ~isempty(number)
        fs_extract = str2double(number{1});
        disp(['提取的采样率是: ', num2str(fs_extract)]);
    else
        disp('未找到符合条件的数字');
        break;
    end
    close all; % 清除图片

    % 遍历文件夹中的每个文件
    for i = 1:length(files)
        % 获取每个文件的完整路径
        fileName = files(i).name;
        fileContent = fileName;

        % 检查文件内容是否以指定的前缀开头
        if strncmp(fileContent, targetPrefix1, length(targetPrefix1))
            disp('匹配成功，执行后续操作');
            % 在此处执行后续代码
            filePath = fullfile(folderPath, fileName);

            % 读取文件内容（根据文件格式调整读取方法）
            fileID = fopen(filePath, 'r');
            fileData = fread(fileID, 'int16')';
            fclose(fileID);

            % 对文件数据进行处理
            disp(['Processing file: ', fileName]);

            idata = fileData(1:2:end-1);            % 复信号
            qdata = fileData(2:2:end);
            sdata = idata+1i*qdata;
            

            % 信号长度
            N = length(sdata);
            % 自动设置窗口大小和重叠样本数
            window = floor(N / 10);  % 设置窗口长度为信号长度的1/10
            noverlap = floor(window / 2);  % 设置重叠样本数为窗口的一半
            nfft = 2^nextpow2(window);  % 设置FFT点数为2的下一个最接近的整数次幂

            % 确保窗口长度不大于信号长度
            if window > N
                window = N;
                noverlap = N - 1;
            end

            % 绘制波形图
            figure;
            subplot(2,1,1);plot(real(sdata));title("实部和虚部波形图");subplot(2,1,2);plot(imag(sdata));
            sdata = sdata.*exp(-1i*2*pi*(fs_extract/2)*(1:length(sdata))/fs_extract);
            % [s, f, t] = spectrogram(sdata,window,noverlap,nfft,fs_extract);
            figure; spectrogram(sdata,window,noverlap,nfft,fs_extract,'yaxis');title('Time-Frequency Plot');
        % else
        %     disp('匹配失败');
        end
    end

    

end