function eta = select_eta(gammat,Channel,phi,P,K,L,sigmat2,sigmar2,varargin)
% select_eta - 选择雷达功率分配系数eta
%   eta = select_eta(gammat,Channel,phi,P,K,L,sigmat2,sigmar2)
%   选择最小的eta值，使得雷达SNR满足要求gammat
%   
%   参数:
%       gammat - 雷达SNR门限值
%       Channel - 信道参数结构体
%       phi - RIS相位矩阵
%       P - 总发射功率
%       K - 用户数量
%       L - 雷达脉冲数
%       sigmat2 - 目标反射系数方差
%       sigmar2 - 雷达噪声功率
%   
%   返回值:
%       eta - 雷达功率分配系数 [0,1]

% 默认不计入通信对雷达的干扰，保持历史结果兼容
includeCommInterferenceInRadar = false;
radarMode = 'nsp';
if ~isempty(varargin)
    v1 = varargin{1};
    if ischar(v1) || isstring(v1)
        radarMode = char(v1);
    else
        includeCommInterferenceInRadar = logical(v1);
    end
end
if numel(varargin) >= 2
    v2 = varargin{2};
    if ischar(v2) || isstring(v2)
        radarMode = char(v2);
    else
        includeCommInterferenceInRadar = logical(v2);
    end
end

% 计算等效用户信道
Hk = Channel.Hu + Channel.Hru*diag(phi)*Channel.G;

% 定义eta的搜索范围和步长（0到1，步长0.05）
etas = 0:0.05:1;

% 初始化为最大值
eta = etas(end);
is_feasible = false;

% 从最小的eta开始搜索，找到满足雷达SNR要求的最小eta
for e = etas
    % 为当前eta设计波束赋形向量
    [Wc,wr] = design_w(Hk,Channel.hdt,Channel.hrt,Channel.G,phi,P,K,e,radarMode);
    
    % 计算雷达SNR
    g = radar_snr(Channel.hdt,Channel.hrt,Channel.G,phi,wr,L,sigmat2,sigmar2, ...
        Wc,includeCommInterferenceInRadar);
    
    % 如果满足要求，更新eta并退出循环
    if g >= gammat
        eta = e;
        is_feasible = true;
        break;
    end
end

if ~is_feasible
    warning('select_eta:InfeasibleConstraint', ...
        '感知约束不可行：eta=1时仍无法满足 gammat=%.4g。已返回 eta=1。', gammat);
end
end