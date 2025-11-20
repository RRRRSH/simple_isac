function phi = compute_phi(mode,Channel)
% compute_phi - 计算RIS相位矩阵
%   phi = compute_phi(mode,Channel) 返回根据指定模式计算的RIS相位矩阵
%   
%   参数:
%       mode - 相位计算模式: 
%             'fixed_target' - 基于目标反射的固定相位模式
%             'random' - 随机相位模式
%       Channel - 信道参数结构体，包含以下字段:
%             .hrt - RIS-目标信道向量 (N×1)
%             .G - RIS-基站信道矩阵 (N×M)
%             .hdt - 基站-目标直接信道向量 (M×1)
%   
%   返回值:
%       phi - RIS相位矩阵 (N×1)，元素为单位模长的复数

% 获取RIS单元数量
N = length(Channel.hrt);

% 根据不同模式选择相位计算方法
switch mode
    case 'fixed_target'
        % 固定目标模式：优化雷达探测性能的相位设计
        % 基于雷达信号路径的相位对准原理
        b = 2*diag(Channel.hrt')*conj(Channel.G)*Channel.hdt;
        phi = exp(1i*angle(b));  % 仅保留相位信息，模长为1
        
    case 'random'
        % 随机模式：生成均匀分布的随机相位
        phi = exp(1i*2*pi*rand(N,1));  % 相位在[0,2π)均匀分布
        
    otherwise
        % 默认情况：所有相位为1（0相位）
        phi = ones(N,1);
end
end