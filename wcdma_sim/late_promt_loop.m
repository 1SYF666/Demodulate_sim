for i = 1:num_symbol
    rr_dec6 = r_dec6( 1 + despreadlen * (i-1) : despreadlen * i );
    rr_dec6_I = real(rr_dec6);
    rr_dec6_Q = imag(rr_dec6);

    Earlyi = circshift(pn_sps_i,-floor(sps/2));
    Earlyq = circshift(pn_sps_q,-floor(sps/2));
    Prompti = circshift(pn_sps_i, 0);
    Promptq = circshift(pn_sps_q, 0);
    Latei = circshift(pn_sps_i,floor(sps/2));
    Lateq = circshift(pn_sps_q,floor(sps/2));

    % 超前
    i_E = rr_dec6_I.*Earlyi;q_E = rr_dec6_Q.*Earlyq;
    I_E = sum(i_E); Q_E = sum(q_E);
    E(i) = abs(sqrt(I_E.^2 + Q_E.^2));
    % 即时
    i_P = rr_dec6_I.*Prompti; q_P = rr_dec6_Q.*Promptq;
    I_P = sum(i_P);Q_P = sum(q_P);
    P(i) = abs(sqrt(I_P.^2));
    output_tra(1 + despreadlen * (i-1) : despreadlen * i) = i_P + 1i * q_P;
    % 滞后
    i_L = rr_dec6_I.*Latei; q_L = rr_dec6_Q.*Lateq;
    I_L = sum(i_L); Q_L = sum(q_L);
    L(i) = abs(sqrt(I_L.^2+Q_L.^2));
    % 环路
    j = i+1;
    loop_error_In(j) = 1/2*(E(i)^2 - L(i)^2)/(E(i)^2 + L(i)^2);
    PLL_Phase_Part(j)=loop_error_In(j)*C11;
    PLL_Freq_Part(j)=loop_error_In(j)*C22+PLL_Freq_Part(j-1);
    loop_error_Out(j)=PLL_Phase_Part(j)+PLL_Freq_Part(j);
    temp2_diff = loop_error_Out(j)-loop_error_Out(j-1);

    % 比较EPL大小
    if L(i)>P(i)&&E(i)<P(i)
        % 信号滞后PN码一个点
        pn_sps_i_temp = circshift(pn_sps_i,floor(sps/4));
        pn_sps_i = pn_sps_i_temp;
        pn_sps_q_temp = circshift(pn_sps_q,floor(sps/4));
        pn_sps_q = pn_sps_q_temp;

    end
    if L(i)<P(i)&&E(i)>P(i)
        % 信号超前PN码一个点
        pn_sps_i_temp = circshift(pn_sps_i,-floor(sps/4));
        pn_sps_i = pn_sps_i_temp;
        pn_sps_q_temp = circshift(pn_sps_q,floor(sps/4));
        pn_sps_q = pn_sps_q_temp;
    end
    
end