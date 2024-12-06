%%
clc;
close all;
clear;
%%
% 初始化 N =2
H_initial = [1 1;1 -1];
N = 4; % 阶数
i = 1;
while i <=log2(N) 
    if i == 1
        H = H_initial;
    else
        H = [H H;H -H];
    end
    i = i+1;
end
 
%% 自相关
% for k = 2 : length(H)/4
%     for ii = 1 :length(H(k,:))*10
%         xcorrout(ii) = sum(H(k,:).*circshift(H(k,:),ii-1))/length(H);
%     end
%     
%     figure;plot(xcorrout);axis([1 length(H(k,:))*3 -0.5 1.2]);grid on;
%     title("自相关");
% end

%% 互相关
for k = 4 : length(H)
    for ii = 1 :length(H(k,:))*10
        xcorrout2(ii) = sum(H(3,:).*circshift(H(k,:),ii-1))/length(H);
    end
    figure;plot(xcorrout2);axis([1 length(H(k,:))*3 -1.2 1.2]);grid on;
    title("互相关");
end


 
 