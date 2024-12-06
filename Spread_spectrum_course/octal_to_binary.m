% 八进制转化为二进制
function [Binary,Binary_temp] = octal_to_binary(number_octal,octal,binary)
    % octal = 8;
    % binary = 2;
    % number_octal = 103;
    a = num2str(number_octal);
    for i = 1 : length(a)
        b(i) = str2num(a(i));
    end
    % 八进制到十进制
    integer = 0;
    N1 = length(b);
    for i = 1 : N1
        integer = integer + b(i)* octal^(N1-i);
    end

    % 十进制到二进制
    j = 0;
    temp = integer;
    while 1
        j = j + 1;
        Binary_temp(j) = mod(temp,binary); % 取余
        temp = ( temp- Binary_temp(j) )/ binary; % 取商
        if temp == 0
            break;
        end
    end
    % 翻转
    N2 = length(Binary_temp);
    for i = 1 : N2
        Binary(i) = Binary_temp(N2-i+1);
    end
end