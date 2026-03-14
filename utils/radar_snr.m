function g = radar_snr(hdt,hrt,G,phi,wr,L,sigmat2,sigmar2,varargin)
% radar_snr - 计算RIS辅助的雷达SNR
%   g = radar_snr(hdt,hrt,G,phi,wr,L,sigmat2,sigmar2)
%   计算考虑RIS反射的雷达信噪比
%   g = radar_snr(...,Wc,true) 计入通信波束对雷达的干扰
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
%       Wc - (可选) 通信波束赋形矩阵 (M×K)
%       includeCommInterference - (可选) 是否计入通信对雷达干扰，默认false
%   
%   返回值:
%       g - 雷达信噪比

% 默认不计入通信对雷达的干扰，保持历史结果兼容
Wc = [];
includeCommInterference = false;
if ~isempty(varargin)
	Wc = varargin{1};
	if numel(varargin) >= 2
		includeCommInterference = varargin{2};
	else
		includeCommInterference = true;
	end
end

h_forward = hdt.' + hrt.'*diag(phi)*G; % 基站 -> 目标 (1xM)
a_target = h_forward * wr;             % 目标处雷达信号（可为1xNr）

% 假设单站雷达，回程信道是去程的转置
h_backward = h_forward.';              % 目标 -> 基站 (Mx1)

% 信号项：多个雷达波束时按总发射功率叠加
signal = L * sigmat2 * sum(abs(a_target).^2) * (norm(h_backward)^2);

% 干扰项：可选计入通信波束回波功率
interference = 0;
if includeCommInterference && ~isempty(Wc)
	interference = L * sigmat2 * sum(abs(h_forward*Wc).^2) * (norm(h_backward)^2);
end

g = real(signal / (sigmar2 + interference));