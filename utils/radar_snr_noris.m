function g = radar_snr_noris(hdt,wr,L,sigmat2,sigmar2,varargin)
% radar_snr_noris - 计算无RIS情况下的雷达SNR
%   g = radar_snr_noris(hdt,wr,L,sigmat2,sigmar2)
%   计算不考虑RIS反射的雷达信噪比
%   g = radar_snr_noris(...,Wc,true) 计入通信波束对雷达的干扰
%   
%   参数:
%       hdt - 基站-目标直接信道 (M×1)
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

% 无RIS时，去程等效信道仅为直接路径
h_forward = hdt.';
a_target = h_forward*wr;

% 单站雷达下，回程信道与去程互易
h_backward = h_forward.';

signal = L*sigmat2*sum(abs(a_target).^2)*(norm(h_backward)^2);
interference = 0;
if includeCommInterference && ~isempty(Wc)
	interference = L*sigmat2*sum(abs(h_forward*Wc).^2)*(norm(h_backward)^2);
end
g = real(signal/(sigmar2 + interference));
end