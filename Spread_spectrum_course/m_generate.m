% m1 进行双极性变换
% Binary_temp：c0 c1 c2 c3 ...
function [m1,m2,len] = m_generate(Binary_temp)
    c = Binary_temp(2:end);               
    N = length(c);                        % c1 ... cN  
    init_register = zeros(1,N);
    init_register(N) = 1;
    temp_register = init_register;
    len = 2^N -1;
    for i = 1:len
        m2(i) = temp_register(N);         % 输出
        and_out = mod(sum(c.*temp_register),2);
        temp_register(2:end) = temp_register(1:end-1);
        temp_register(1) = and_out;
    end
    m1 = 1 - m2*2; % 1-->-1；0-->1
end