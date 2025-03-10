%%%%% 参数配置
function out = GMSK_Mod(inBits,OSR,Tb,BT)


% inBits = zeros(1,len);
in = 1-2*inBits;
[out_I,out_Q] = gmsk_mod(in,Tb,OSR,BT);
out = out_I+1j*out_Q;
end


function [I,Q] = gmsk_mod(BURST,Tb,OSR,BT)
%
% gmsk_mod:   This function accepts a GSM burst bit sequence and
%             performs a GMSK modulation of the sequence. The
%             modulation is according to the GSM 05.05 recommendations
%
% SYNTAX:     [i,q] = gmsk_mod(burst,Tb,osr,BT)
%
% INPUT:      burst   A differential encoded bit sequence (-1,+1)
%             Tb      Bit duration (GSM: Tb = 3.692e-6 Sec.)
%             osr     Simulation oversample ratio. osr determines the
%                     number of simulation steps per information bit
%             BT      The bandwidth/bit duration product (GSM: BT = 0.3)
%
% OUTPUT:     i,q     In-phase (i) and quadrature-phase (q) baseband
%                     representation of the GMSK modulated input burst 
%                     sequence
%
% SUB_FUNC:   ph_g.m  This sub-function is required in generating the
%                     frequency and phase pulse functions.
%
% WARNINGS:   Sub-function ph_g.m assumes a 3xTb frequency pulse 
%             truncation time
%
% TEST(S):    Function tested using the following relations
%             
%             i)
%
%             I^2 + Q^2 = Cos(a)^2 + Sin(a)^2 = 1
%
%             ii)
%
%             When the input consists of all 1's the resulting baseband
% 	      outputs the function should return a sinusoidal signal of
%             frequency rb/4, i.e. a signal having a periode time of 
%             approximately 4*Tb = 4*3.692e-6 s = 1.48e-5 s for GSM
%
% AUTHOR:   Jan H. Mikkelsen / Arne Norre Ekstr鴐
% EMAIL:    hmi@kom.auc.dk / aneks@kom.auc.dk
%
% $Id: gmsk_mod.m,v 1.5 1998/02/12 10:50:10 aneks Exp $

% ACCUIRE GMSK FREQUENCY PULSE AND PHASE FUNCTION
%
[g,~] = ph_g(Tb,OSR,BT);

% PREPARE VECTOR FOR DATA PROCESSING
%
bits = length(BURST);           % 长度为148bit
f_res = zeros(1,(bits+4)*OSR);

% GENERATE RESULTING FREQUENCY PULSE SEQUENCE
%
for n = 1:bits
%     f_res((n-2)*OSR+1:(n+3)*OSR) = f_res((n-2)*OSR+1:(n+3)*OSR) + BURST(n).*g;
  f_res((n-1)*OSR+1:(n+4)*OSR) = f_res((n-1)*OSR+1:(n+4)*OSR) + BURST(n).*g;
end

% CALCULATE RESULTING PHASE FUNCTION
%
theta = pi*cumsum(f_res);

% PREPARE DATA FOR OUTPUT
%
I = cos(theta);
Q = sin(theta);

end


function [G_FUN, Q_FUN] = ph_g(Tb,OSR,BT)
%
% PH_G:     This function calculates the frequency and phase functions
%           required for the GMSK modulation. The functions are 
%           generated according to the GSM 05.05 recommendations
%
% SYNTAX:   [g_fun, q_fun] = ph_g(Tb,osr,BT)
%
% INPUT:    Tb      Bit duration (GSM: Tb = 3.692e-6 Sec.)
%           osr     Simulation oversample ratio. osr determines the
%                   number of simulation steps per information bit
%           BT      The bandwidth/bit duration product (GSM: BT = 0.3)
%
% OUTPUT:   g_fun, q_fun   Vectors contaning frequency and phase 
%                          function outputs when evaluated at osr*tb
%
% SUB_FUNC: None
%
% WARNINGS: Modulation length of 3 is assumed !
%
% TEST(S):  Tested through function gsmk_mod.m
%
% AUTHOR:   Jan H. Mikkelsen / Arne Norre Ekstr鴐
% EMAIL:    hmi@kom.auc.dk / aneks@kom.auc.dk
%
% $Id: ph_g.m,v 1.6 1998/02/12 10:50:54 aneks Exp $

% SIMULATION SAMPLE FREQUENCY
%
Ts = Tb/OSR;

% PREPARING VECTORS FOR DATA PROCESSING
%
PTV = -2*Tb:Ts:2*Tb;
RTV = -Tb/2:Ts:Tb/2-Ts;

% GENERATE GAUSSIAN SHAPED PULSE
%
sigma = sqrt(log(2))/(2*pi*BT);
gauss = (1/(sqrt(2*pi)*sigma*Tb))*exp(-PTV.^2/(2*sigma^2*Tb^2));   % 写出高斯函数

% GENERATE RECTANGULAR PULSE
%
rect = 1/(2*Tb)*ones(size(RTV));

% CALCULATE RESULTING FREQUENCY PULSE
%
G_TEMP = conv(gauss,rect);

% TRUNCATING THE FUNCTION TO 3xTb
%
% G = G_TEMP(OSR+1:4*OSR); 
G = G_TEMP;

% TRUNCATION IMPLIES THAT INTEGRATING THE FREQUENCY PULSE
% FUNCTION WILL NOT EQUAL 0.5, HENCE THE RE-NORMALIZATION
%
G_FUN = (G-G(1))./(2*sum(G-G(1)));

% CALCULATE RESULTING PHASE PULSE
%
Q_FUN = cumsum(G_FUN);
end


