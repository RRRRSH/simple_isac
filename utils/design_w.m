function [Wc,wr] = design_w(Hk,hdt,P,K,eta)
% design_w - 设计通信和雷达波束赋形向量
%   [Wc,wr] = design_w(Hk,hdt,P,K,eta) 计算通信和雷达的波束赋形向量
%   
%   参数:
%       Hk - 等效用户信道矩阵 (K×M)
%       hdt - 基站-目标直接信道向量 (M×1)
%       P - 总发射功率
%       K - 用户数量
%       eta - 雷达功率分配系数 [0,1]，表示分配给雷达的功率比例
%   
%   返回值:
%       Wc - 通信波束赋形矩阵 (M×K)
%       wr - 雷达波束赋形向量 (M×1)

% 使用Moore-Penrose伪逆设计通信波束赋形矩阵（ZF预编码）
Wc = Hk'/(Hk*Hk');

% 对每个用户的波束赋形向量进行归一化
for k = 1:K
    Wc(:,k) = Wc(:,k)/norm(Wc(:,k)+1e-12);  % 添加小量避免除零
end

% 计算分配给通信的功率
Pc = (1-eta)*P;  % 通信功率 = 总功率 * (1-eta)

% 分配功率给通信波束赋形矩阵
Wc = Wc*diag(sqrt(Pc/K)*ones(1,K));  % 均匀分配功率给每个用户

% 设计雷达波束赋形向量（最大比传输）
wr = hdt/norm(hdt+1e-12)*sqrt(eta*P);  % 归一化后分配雷达功率
end