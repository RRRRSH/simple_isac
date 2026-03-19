function [phi_joint,Wc_joint,Wr_joint,Vsr_joint,gammat_joint] = run_joint_snr_optimization(Channel,P,K,L,sigmat2,sigmar2,sigmak2,gammat,varargin)
% run_joint_snr_optimization - 使用SNR约束联合优化同时更新RIS相位与发射波束
%   输出通信波束Wc_joint与雷达波束矩阵Wr_joint（可为多列）

M = size(Channel.Hu,2);
N = length(Channel.hrt);

if nargin >= 10 && ~isempty(varargin{1}) && ~isempty(varargin{2})
    phi0 = varargin{1};
    W0 = varargin{2};
    Wc0 = W0(:, 1:K);
    wr0 = W0(:, K+1:end);
    Hk0 = Channel.Hu + Channel.Hru*diag(phi0)*Channel.G;
else
    % 先运行一次完全严格的波束凸优化，得到针对当前 phi0 的理论最优稳健波束
    % 防止只给一个 eta=0.5 的初值由于约束收得太紧导致 update phi 第一轮时无解成 NaN
    phi0 = compute_phi('fixed_target',Channel);
    Hk0 = Channel.Hu + Channel.Hru*diag(phi0)*Channel.G;
    [Wc0,wr0] = optimize_w_for_fixed_phi(Channel,phi0,P,K,L,sigmat2,sigmar2,sigmak2,gammat);

    % wr0 已经是包含雷达维度 (M x M) 的完整矩阵，不需要再像单列向量模式那样补零
    W0 = [Wc0, wr0];
end

Prms.M = M;
Prms.N = N;
Prms.K = K;
Prms.sigmar2 = sigmar2;
Prms.sigmak2 = sigmak2;
Prms.sigmat2 = sigmat2;
Prms.Nmax = 100;
Prms.res_th = 5e-4;
Prms.L = L;
Prms.gammat = gammat;
Prms.P = P;

try
    [W_joint,phi_joint,Vsr_joint,gammat_joint] = get_W_phi_SNR(Prms,Channel,phi0,W0);
catch ME
    cvx_clear;
    fprintf('\n=> [FATAL ERROR] 联合优化抛出异常: %s\n', ME.message);
    fprintf('=>        [Stack] File: %s, Line: %d\n\n', ME.stack(1).file, ME.stack(1).line);
    W_joint = W0;
    phi_joint = phi0;
    Vsr_joint = sum_rate(Hk0,Wc0,sigmak2,wr0,false);
    gammat_joint = radar_snr(Channel.hdt,Channel.hrt,Channel.G,phi0,wr0,L,sigmat2,sigmar2,Wc0,false);
end

Wc_joint = W_joint(:,1:K);
% 论文联合优化中所有流均可用于感知，因此雷达评估应使用完整W矩阵
Wr_joint = W_joint;
end
