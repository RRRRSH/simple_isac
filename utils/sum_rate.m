function R = sum_rate(Hk,Wc,sigmak2)
% sum_rate - 计算多用户通信系统的总和速率
%   R = sum_rate(Hk,Wc,sigmak2) 计算所有用户的总频谱效率
%   
%   参数:
%       Hk - 用户信道矩阵 (K×M)
%       Wc - 通信波束赋形矩阵 (M×K)
%       sigmak2 - 用户噪声功率
%   
%   返回值:
%       R - 总和速率 (bps/Hz)

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
    
    % 计算用户k的速率（香农容量）并累加
    R = R + log2(1 + num/den);
end
end