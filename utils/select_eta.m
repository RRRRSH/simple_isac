function eta = select_eta(gammat,Channel,phi,P,K,L,sigmat2,sigmar2)
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

% 计算等效用户信道
Hk = Channel.Hu + Channel.Hru*diag(phi)*Channel.G;

% 定义eta的搜索范围和步长（0到0.5，步长0.05）
etas = 0:0.05:0.5;

% 初始化为最大值
eta = etas(end);

% 从最小的eta开始搜索，找到满足雷达SNR要求的最小eta
for e = etas
    % 为当前eta设计波束赋形向量
    [~,wr] = design_w(Hk,Channel.hdt,P,K,e);
    
    % 计算雷达SNR
    g = radar_snr(Channel.hdt,Channel.hrt,Channel.G,phi,wr,L,sigmat2,sigmar2);
    
    % 如果满足要求，更新eta并退出循环
    if g >= gammat
        eta = e;
        break;
    end
end
end