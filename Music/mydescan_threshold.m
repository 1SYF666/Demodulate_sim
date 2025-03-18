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
data = diag(Lambda);
[row,col] = size(data);
% 运行DBSCAN
eps    = 0.5;    % 邻域半径
minPts = 5;      % 最小邻域点数
labels = myDBSCAN(data, eps, minPts);

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

% 绘制噪音点(如果有)
idxNoise = (labels==0);

if col == 1
    scatter(data(idxNoise,1), 0, 36, 'k', 'x');
elseif col == 2
    scatter(data(idxNoise,1), data(idxNoise,2), 36, 'k', 'x');
end

title('DBSCAN Clustering Result');
xlabel('X'); ylabel('Y');
