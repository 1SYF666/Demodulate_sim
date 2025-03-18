function labels = myDBSCAN(data, eps, minPts)
%MYDBSCAN  实现基于密度的聚类算法(DBSCAN)
%
%   labels = myDBSCAN(data, eps, minPts)
%
%   输入:
%       data   - m x d 的数据矩阵, m 为样本数, d 为特征维度
%       eps    - 距离阈值 (ε-邻域半径)
%       minPts - 形成核心点的最少邻域点数
%
%   输出:
%       labels - size为 m x 1 的聚类标签向量
%                - 未分到任何簇的点标记为 0
%                - 第 k 个簇标记为整数 k
%
%   参考步骤:
%   1) 从数据集中任意选取一个未处理的点 p
%   2) 若 p 为核心点 (邻域内点数 >= minPts), 则以 p 为种子扩展出一个簇
%   3) 若 p 是边界点 (邻域内点数 < minPts), 则将其标记为噪音或暂留
%   4) 重复, 直至所有点都被处理过

[m, ~] = size(data);     % m 为样本数
labels = zeros(m, 1);    % 初始化所有点的聚类标签为 0 (表示未分配 / 噪音)
clusterId = 0;           % 用于记录当前正在创建的簇编号

% 用于标记点的访问状态, 以免重复处理. 0=未访问, 1=已访问
visited = false(m, 1);

for i = 1 : m
    % 若该点已访问, 则跳过
    if visited(i)
        continue;
    end
    
    % 否则, 标记为已访问
    visited(i) = true;
    
    % 查找点 i 在 eps 范围内的邻域点索引
    neighbors = regionQuery(data, i, eps);
    
    % 判断点 i 是否为核心点
    if numel(neighbors) < minPts
        % 邻域太小, i 被标记为噪音 (或可能是边界点)
        labels(i) = 0; % 可以暂时记为 0
    else
        % 找到了一个核心点, 创建一个新簇
        clusterId = clusterId + 1;
        % 将 i 点及其邻域都扩展到该簇
        labels = expandCluster(data, labels, i, neighbors, clusterId, eps, minPts, visited);
    end
end

end % end of myDBSCAN


%% ======================== 子函数 1: expandCluster =======================
function labels = expandCluster(data, labels, seedIndex, neighbors, clusterId, eps, minPts, visited)
% 将核心点 seedIndex 的邻域加入聚类 clusterId, 并继续向外扩展

% 先把 seedIndex 这个核心点分配到当前簇
labels(seedIndex) = clusterId;

% 用一个队列来处理待扩展的点. 初始时包含所有邻域点
searchQueue = neighbors;

idx = 1; % 指向 searchQueue 中当前要处理的位置
while idx <= numel(searchQueue)
    % 取出一个待处理点
    currentPoint = searchQueue(idx);
    
    % 如果该点还没访问过, 则标记访问
    if ~visited(currentPoint)
        visited(currentPoint) = true;
        
        % 获取该点的邻域
        currentNeighbors = regionQuery(data, currentPoint, eps);
        
        % 如果该点也是一个核心点(邻域>=minPts),
        % 则把它的邻域也加入搜索队列, 以便继续向外扩张
        if numel(currentNeighbors) >= minPts
            searchQueue = [searchQueue; currentNeighbors]; %#ok<AGROW>
        end
    end
    
    % 如果该点还未被分配到某个簇 (标签=0), 分配到当前簇
    if labels(currentPoint) == 0
        labels(currentPoint) = clusterId;
    end
    
    idx = idx + 1;
end

end % end of expandCluster


%% ======================== 子函数 2: regionQuery =========================
function neighbors = regionQuery(data, idx, eps)
% 返回与 data(idx,:) 的欧氏距离 <= eps 的所有点索引

distVec = sqrt(sum((data - data(idx,:)) .^ 2, 2)); 
neighbors = find(distVec <= eps);

end
