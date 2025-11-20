function g = radar_snr(hdt,hrt,G,phi,wr,L,sigmat2,sigmar2)
% radar_snr - 计算RIS辅助的雷达SNR
%   g = radar_snr(hdt,hrt,G,phi,wr,L,sigmat2,sigmar2)
%   计算考虑RIS反射的雷达信噪比
%   
%   参数:
%       hdt - 基站-目标直接信道 (M×1)
%       hrt - RIS-目标信道 (N×1)
%       G - RIS-基站信道 (N×M)
%       phi - RIS相位矩阵 (N×1)
%       wr - 雷达波束赋形向量 (M×1)
%       L - 雷达脉冲数
%       sigmat2 - 目标反射系数方差
%       sigmar2 - 雷达噪声功率
%   
%   返回值:
%       g - 雷达信噪比

% 计算组合信道响应（直接路径 + RIS反射路径）
a = (hdt.' + hrt.'*diag(phi)*G)*wr;

% 计算雷达SNR（考虑相干积累L个脉冲）
g = real(L*sigmat2*abs(a)^2/sigmar2);
end