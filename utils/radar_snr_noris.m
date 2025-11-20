function g = radar_snr_noris(hdt,wr,L,sigmat2,sigmar2)
% radar_snr_noris - 计算无RIS情况下的雷达SNR
%   g = radar_snr_noris(hdt,wr,L,sigmat2,sigmar2)
%   计算不考虑RIS反射的雷达信噪比
%   
%   参数:
%       hdt - 基站-目标直接信道 (M×1)
%       wr - 雷达波束赋形向量 (M×1)
%       L - 雷达脉冲数
%       sigmat2 - 目标反射系数方差
%       sigmar2 - 雷达噪声功率
%   
%   返回值:
%       g - 雷达信噪比

% 计算直接路径信道响应
% 无RIS时，只有基站到目标的直接路径
a = hdt.'*wr;

% 计算雷达SNR
g = real(L*sigmat2*abs(a)^2/sigmar2);  
end