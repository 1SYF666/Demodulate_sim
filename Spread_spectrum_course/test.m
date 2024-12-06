clc;close all;clear;
%% m序列生成


integer_Octal1 = 147;
[binary_temp1] = m_generate(integer_Octal1);
ff = binary_temp1(2:end); 

integer_Octal2 = 155;
[binary_temp2] = m_generate(integer_Octal2);
ff2 = binary_temp2(2:end); 

N = length(ff);
init_register = zeros(1,N);
init_register(N) = 1;
temp_register = init_register;

init_register2 = zeros(1,N);
init_register2(N-1) = 1;
temp_register2 = init_register2;

length_m = 2^N -1;
for i = 1:length_m
   msequence(i) = temp_register(N);       % 序列一输出
   and_out = mod(sum(ff.*temp_register),2);
   temp_register(2:end) = temp_register(1:end-1);
   temp_register(1) = and_out; 
   
   msequence2(i) = temp_register2(N);       % 序列二输出
   and_out2 = mod(sum(ff2.*temp_register2),2);
   temp_register2(2:end) = temp_register2(1:end-1);
   temp_register2(1) = and_out2; 
   
   gold_msequence(i) = mod(msequence(i)+msequence2(i),2);
end

% 自相关函数

m_mesquence1 = 1 - gold_msequence*2;

for ii = 1 :length(m_mesquence1)*10
    xcorrout(ii) = sum(m_mesquence1.*circshift(m_mesquence1,ii-1))/length_m;
end

figure;plot(xcorrout);axis([1 length(m_mesquence1)*3 -0.5 2]);grid on;
title("自相关");


%% 另一个gold
integer_Octal1 = 155;
[binary_temp1] = m_generate(integer_Octal1);
ff = binary_temp1(2:end); 

integer_Octal2 = 103;
[binary_temp2] = m_generate(integer_Octal2);
ff2 = binary_temp2(2:end); 

N = length(ff);
init_register = zeros(1,N);
init_register(N) = 1;
temp_register = init_register;

init_register2 = zeros(1,N);
init_register2(N-1) = 1;
temp_register2 = init_register2;

length_m = 2^N -1;
for i = 1:length_m
   msequence(i) = temp_register(N);       % 序列一输出
   and_out = mod(sum(ff.*temp_register),2);
   temp_register(2:end) = temp_register(1:end-1);
   temp_register(1) = and_out; 
   
   msequence2(i) = temp_register2(N);       % 序列二输出
   and_out2 = mod(sum(ff2.*temp_register2),2);
   temp_register2(2:end) = temp_register2(1:end-1);
   temp_register2(1) = and_out2; 
   
   gold_msequence(i) = mod(msequence(i)+msequence2(i),2);
end

% 自相关函数

m_mesquence2 = 1 - gold_msequence*2;

for ii = 1 :length(m_mesquence2)*10
    xcorrout(ii) = sum(m_mesquence2.*circshift(m_mesquence2,ii-1))/length_m;
end

figure;plot(xcorrout);axis([1 length(m_mesquence2)*3 -0.5 2]);grid on;
title("自相关");


%% 互相关

for ii = 1 :length(m_mesquence2)*10 
    xcorrout3(ii) = sum(m_mesquence1.*circshift(m_mesquence2,ii-1))/length_m;
end

figure;plot(xcorrout3);

%% 函数
function [binary_temp] = m_generate(integer_Octal)
    a = num2str(integer_Octal); 
    for i = 1 : length(a)
        b(i) = str2num(a(i));
    end

    % 整数_to_二进制 ;
    for i = 1 : length(b)
        j = 0;
        temp = b(i);
        while 1
            j = j + 1;
            binary(i,j) = mod(temp,2); % 取余
            temp = ( temp- binary(i,j) )/ 2; % 取商
            if temp == 1
                binary(i,j+1) = 1;
                break;
            end
            if temp == 0
                break;
            end
        end
    end

    % 翻转
    [row,col] = size(binary);

    k = 0;
    for i = 1 : 1 : row
        flag1 = 1;
        flag2 = 0;
        for j = col : -1 : 1 
            if flag1
                if binary(i,j)
                    flag1 = 0;
                    flag2 = 1;
                end
            end

            if flag2
                k = k + 1;
                binary_temp(k) = binary(i,j);
            end

        end
    end

end












