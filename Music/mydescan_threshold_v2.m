%% ********* music simulation ********* %%
%% ***** data:20250316 authoor:ShenYifu ****  %%
%{
该程序作用：
    对噪声数据选择聚类最优阈值 eps and minPts
参考资料：
    基于特征值聚类的MUSIC算法
%}
clc;close all;clear;

% 生成一些测试数据
%{
rng(1); % 保持随机一致性
data1 = randn(50,2)*0.2 + [0,0];
data2 = randn(50,2)*0.2 + [2,2];
data3 = randn(50,2)*0.2 + [5,0];
data2D = [data1; data2; data3];
%}

load Lambda.mat;
data = diag(Lambda); % m*d -- m个簇

[row,col] = size(data);
if col == 1
    Reps = [0.2: 0.25 : (max(data(:))- min(data(:)))/2];
    Minpts = [3:1:10];
else
    fprintf("此时输入数据不是N*1,需要仔细考虑Reps的选取！")
    return;
end
for i = 1 : length(Reps)
    for j =  1 : length(Minpts)
        % 运行DBSCAN
        eps    = Reps(i);    % 邻域半径
        minPts = Minpts(j);      % 最小邻域点数
        labels = myDBSCAN(data, eps, minPts);
        idxNoise = (labels==0);% 绘制噪音点(如果有)
        numOnes = sum(idxNoise);         % 逻辑向量为 true (1) 的数量
        numZeros = sum(~idxNoise);       % 逻辑向量为 false (0) 的数量
        if  numOnes > numZeros
            fprintf("(%f,%f)噪声点太多\n",eps,minPts);
            continue;
        elseif numOnes == 0
            fprintf("(%f,%f)没有噪声点\n",eps,minPts);
            continue;
        end
        % 可视化结果
        figure; hold on;
        clusterIds = unique(labels);
        clusterIds(clusterIds==0) = []; % 0 表示噪音点, 单独处理
        colors = lines(length(clusterIds)); % 色表

        for k = 1:numel(clusterIds)
            idxk = (labels == clusterIds(k));
            if col == 1
                scatter(data(idxk,1), 0, 36, colors(k,:), 'filled');
            elseif col == 2
                scatter(data(idxk,1), data(idxk,2), 36, colors(k,:), 'filled');
            end
        end

        if col == 1
            scatter(data(idxNoise,1), 0, 36, 'k', 'x');
        elseif col == 2
            scatter(data(idxNoise,1), data(idxNoise,2), 36, 'k', 'x');
        end
        title(['DBSCAN Clustering of MUSIC Peaks (Reps: ', num2str(eps), ', Minpts: ', num2str(minPts), ')']);
        xlabel('X'); ylabel('Y');

    end
end

