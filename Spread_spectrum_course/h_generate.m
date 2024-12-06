% hadamard矩阵产生
function [H] = h_generate(N)
    H_initial = [1 1;1 -1];
    % N = 4; % 阶数
    i = 1;
    while i <=log2(N)
        if i == 1
            H = H_initial;
        else
            H = [H H;H -H];
        end
        i = i+1;
    end
end