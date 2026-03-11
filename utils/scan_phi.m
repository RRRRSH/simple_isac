function phi = scan_phi(objective,Channel,P,K,L,sigmat2,sigmar2,gammat,sigmak2,varargin)
% scan_phi - 通过扫描优化RIS相位
%   phi = scan_phi(objective,Channel,P,K,L,sigmat2,sigmar2,gammat,sigmak2)
%   通过逐一优化每个RIS单元的相位来最大化目标函数
%   
%   参数:
%       objective - 优化目标: 'sumrate' (通信速率) 或 'snr' (雷达SNR)
%       Channel - 信道参数结构体
%       P - 总发射功率
%       K - 用户数量
%       L - 雷达脉冲数
%       sigmat2 - 目标反射系数方差
%       sigmar2 - 雷达噪声功率
%       gammat - 雷达SNR门限值
%       sigmak2 - 用户噪声功率
%   
%   返回值:
%       phi - 优化后的RIS相位矩阵

% 初始化相位为基于目标的固定相位
phi = compute_phi('fixed_target',Channel);
N = length(phi);  % RIS单元数量
includeRadarInterferenceInRate = false;  % 设为true可计入雷达对通信的干扰
includeCommInterferenceInRadar = false;  % 设为true可计入通信对雷达的干扰
if ~isempty(varargin)
    includeRadarInterferenceInRate = varargin{1};
    if numel(varargin) >= 2
        includeCommInterferenceInRadar = varargin{2};
    end
end

% 相位扫描网格（去除与0等价的2π端点，避免重复评估）
grid = linspace(0,2*pi,13);
grid(end) = [];

% 逐一优化每个RIS单元的相位
for n = 1:N
    best = phi(n);  % 当前单元的最佳相位
    bestval = -inf;  % 最佳目标函数值
    
    % 扫描当前单元的所有可能相位值
    for a = 1:numel(grid)
        % 创建候选相位向量，仅修改当前单元的相位
        cand = phi;
        cand(n) = exp(1i*grid(a));
        
        % 为候选相位选择功率分配系数
        eta = select_eta(gammat,Channel,cand,P,K,L,sigmat2,sigmar2,includeCommInterferenceInRadar);
        
        % 计算等效信道
        Hk = Channel.Hu + Channel.Hru*diag(cand)*Channel.G;
        
        % 设计波束赋形向量
        [Wc,wr] = design_w(Hk,Channel.hdt,Channel.hrt,Channel.G,cand,P,K,eta);
        
        % 根据目标计算相应的性能指标
        if strcmp(objective,'sumrate')
            val = sum_rate(Hk,Wc,sigmak2,wr,includeRadarInterferenceInRate);  % 计算总和速率
        else
            val = radar_snr(Channel.hdt,Channel.hrt,Channel.G,cand,wr,L,sigmat2,sigmar2, ...
                Wc,includeCommInterferenceInRadar);  % 计算雷达SNR
        end
        
        % 更新最佳相位和目标值
        if val > bestval
            bestval = val;
            best = cand(n);
        end
    end
    
    % 更新当前单元的相位为最佳值
    phi(n) = best;
end
end