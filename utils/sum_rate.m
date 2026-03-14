function R = sum_rate(Hk,Wc,sigmak2,varargin)
% sum_rate - 计算多用户通信系统的总和速率
%   R = sum_rate(Hk,Wc,sigmak2) 计算所有用户的总频谱效率（默认不计雷达干扰）
%   R = sum_rate(Hk,Wc,sigmak2,wr,true) 计入雷达波束对通信用户的干扰
%   
%   参数:
%       Hk - 用户信道矩阵 (K×M)
%       Wc - 通信波束赋形矩阵 (M×K)
%       sigmak2 - 用户噪声功率
%       wr - (可选) 雷达波束赋形向量 (M×1)
%       includeRadarInterference - (可选) 是否计入雷达干扰，默认false
%   
%   返回值:
%       R - 总和速率 (bps/Hz)

% 默认不计入雷达干扰，保持与历史结果兼容
wr = [];
includeRadarInterference = false;
if ~isempty(varargin)
    wr = varargin{1};
    if numel(varargin) >= 2
        includeRadarInterference = varargin{2};
    else
        includeRadarInterference = true;
    end
end

% 获取用户数量
K = size(Hk,1);

% 初始化总和速率
R = 0;

% 计算每个用户的速率并累加
for k = 1:K
    % 计算信号功率（期望用户信号）
    num = abs(Hk(k,:)*Wc(:,k))^2;
    
    % 计算干扰加噪声功率
    den = 0;
    for j = 1:K
        if j~=k
            % 累加其他用户的干扰功率
            den = den + abs(Hk(k,:)*Wc(:,j))^2;
        end
    end
    % 添加噪声功率
    den = den + sigmak2;

    % 可选：计入雷达波束泄漏对通信用户的干扰
    if includeRadarInterference && ~isempty(wr)
        % 兼容单雷达波束(Mx1)与多雷达波束(MxNr)两种输入
        den = den + sum(abs(Hk(k,:)*wr).^2);
    end
    
    % 计算用户k的速率（香农容量）并累加
    R = R + log2(1 + num/den);
end
end