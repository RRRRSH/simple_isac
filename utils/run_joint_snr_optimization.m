function [phi_joint,Wc_joint,Wr_joint,Vsr_joint,gammat_joint] = run_joint_snr_optimization(Channel,P,K,L,sigmat2,sigmar2,sigmak2,gammat)
% run_joint_snr_optimization - 使用SNR约束联合优化同时更新RIS相位与发射波束
%   输出通信波束Wc_joint与雷达波束矩阵Wr_joint（可为多列）

M = size(Channel.Hu,2);
N = length(Channel.hrt);

% 用固定相位+启发式波束作为可行初值
phi0 = compute_phi('fixed_target',Channel);
eta0 = select_eta(gammat,Channel,phi0,P,K,L,sigmat2,sigmar2,false);
Hk0 = Channel.Hu + Channel.Hru*diag(phi0)*Channel.G;
[Wc0,wr0] = design_w(Hk0,Channel.hdt,Channel.hrt,Channel.G,phi0,P,K,eta0);

% 论文模型采用K+M个数据流：前K列用于通信，后M列用于雷达
W0 = [Wc0, wr0, zeros(M,max(M-1,0))];

Prms.M = M;
Prms.N = N;
Prms.K = K;
Prms.sigmar2 = sigmar2;
Prms.sigmak2 = sigmak2;
Prms.sigmat2 = sigmat2;
Prms.Nmax = 40;
Prms.res_th = 1e-3;
Prms.L = L;
Prms.gammat = gammat;
Prms.P = P;

try
    [W_joint,phi_joint,Vsr_joint,gammat_joint] = get_W_phi_SNR(Prms,Channel,phi0,W0);
catch ME
    warning('run_joint_snr_optimization:SolverFailed', ...
        '联合优化失败，退化为初值方案。原因: %s', ME.message);
    W_joint = W0;
    phi_joint = phi0;
    Vsr_joint = sum_rate(Hk0,Wc0,sigmak2,wr0,false);
    gammat_joint = radar_snr(Channel.hdt,Channel.hrt,Channel.G,phi0,wr0,L,sigmat2,sigmar2,Wc0,false);
end

Wc_joint = W_joint(:,1:K);
% 论文联合优化中所有流均可用于感知，因此雷达评估应使用完整W矩阵
Wr_joint = W_joint;
end
