function [Wc,wr] = design_w(Hk,hdt,hrt,G,phi,P,K,eta,varargin)
% design_w - 设计通信和雷达波束赋形向量
%   [Wc,wr] = design_w(Hk,hdt,hrt,G,phi,P,K,eta) 计算通信和雷达的波束赋形向量
%   [Wc,wr] = design_w(...,'mrt') 使用传统MRT雷达波束
%   [Wc,wr] = design_w(...,'nsp') 使用NSP(零空间投影)雷达波束（默认）
%   
%   参数:
%       Hk - 等效用户信道矩阵 (K×M)
%       hdt - 基站-目标直接信道向量 (M×1)
%       hrt - RIS-目标信道向量 (N×1)
%       G - RIS-基站信道矩阵 (N×M)
%       phi - RIS相位向量 (N×1)
%       P - 总发射功率
%       K - 用户数量
%       eta - 雷达功率分配系数 [0,1]，表示分配给雷达的功率比例
%   
%   返回值:
%       Wc - 通信波束赋形矩阵 (M×K)
%       wr - 雷达波束赋形向量 (M×1)

radarMode = 'nsp';
if ~isempty(varargin)
    radarMode = lower(string(varargin{1}));
end

% 使用Moore-Penrose伪逆设计通信波束赋形矩阵（ZF预编码）
Wc = Hk'/(Hk*Hk');

% 对每个用户的波束赋形向量进行归一化
for k = 1:K
    Wc(:,k) = Wc(:,k)/(norm(Wc(:,k))+1e-12);  % 添加小量避免除零
end

% 计算分配给通信的功率
Pc = (1-eta)*P;  % 通信功率 = 总功率 * (1-eta)

% 分配功率给通信波束赋形矩阵
Wc = Wc*diag(sqrt(Pc/K)*ones(1,K));  % 均匀分配功率给每个用户

% 基于直达+RIS反射的复合感知信道设计雷达波束
h_eff_t = hdt.' + hrt.'*diag(phi)*G;
a_target = conj(h_eff_t.');
wr_mrt = a_target/(norm(a_target)+1e-12);

switch radarMode
    case "mrt"
        wr_dir = wr_mrt;
    case "nsp"
        M = size(Hk,2);
        Ku = size(Hk,1);
        if M > Ku
            P_perp = eye(M) - Hk' * pinv(Hk*Hk') * Hk;
            wr_raw = P_perp * a_target;
            if norm(wr_raw) > 1e-10
                wr_dir = wr_raw / norm(wr_raw);
            else
                warning('design_w:NSPZeroVector', ...
                    'NSP投影后雷达波束近似为0，回退到MRT。');
                wr_dir = wr_mrt;
            end
        else
            warning('design_w:NSPDofInsufficient', ...
                'NSP需要M>K。当前M=%d, K=%d，回退到MRT。', M, Ku);
            wr_dir = wr_mrt;
        end
    otherwise
        warning('design_w:UnknownRadarMode', ...
            '未知雷达波束模式: %s，回退到NSP。', char(radarMode));
        M = size(Hk,2);
        P_perp = eye(M) - Hk' * pinv(Hk*Hk') * Hk;
        wr_raw = P_perp * a_target;
        if norm(wr_raw) > 1e-10
            wr_dir = wr_raw / norm(wr_raw);
        else
            wr_dir = wr_mrt;
        end
end

wr = wr_dir * sqrt(eta*P);
end