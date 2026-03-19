function [Wc_opt,wr_opt,phi_out,Vsr_opt,gammat_opt] = optimize_w_for_fixed_phi(Channel,phi_fixed,P,K,L,sigmat2,sigmar2,sigmak2,gammat,varargin)
% optimize_w_for_fixed_phi
% 在给定RIS相位phi_fixed下，按SNR约束最大化通信总速率：
%   max_{Wc,wr} sum_k log2(1+SINR_k)
%   s.t. SNR_rad >= gammat,
%        ||Wc||_F^2 + ||wr||_2^2 <= P

M = size(Channel.Hu,2);
N = length(Channel.hrt);
phi_out = phi_fixed;
W0_in = [];
verbose = false;
if numel(varargin) >= 1 && ~isempty(varargin{1})
    W0_in = varargin{1};
end
if numel(varargin) >= 2 && ~isempty(varargin{2})
    verbose = logical(varargin{2});
end

if isempty(W0_in)
    % 使用均匀功率分配的启发式解作为初值，提供稳定的CVX启动点
    eta0 = 0.5; % 固定分配50%功率给雷达
    Hk0 = Channel.Hu + Channel.Hru*diag(phi_fixed)*Channel.G;
    [Wc0,wr0] = design_w(Hk0,Channel.hdt,Channel.hrt,Channel.G,phi_fixed,P,K,eta0);
    W0 = [Wc0, wr0, zeros(M,max(M-1,0))];
else
    W0 = W0_in;
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
Prms.freezePhi = true;
Prms.verbose = verbose;

try
    [W_opt,~,Vsr_opt,gammat_opt] = get_W_phi_SNR(Prms,Channel,phi_fixed,W0);
catch ME
    cvx_clear;
    warning('optimize_w_for_fixed_phi:SolverFailed', ...
        '固定phi下的P1波束优化失败，回退启发式初值。原因: %s', ME.message);
    W_opt = W0;
    Vsr_opt = sum_rate(Hk0,Wc0,sigmak2,wr0,false);
    gammat_opt = radar_snr(Channel.hdt,Channel.hrt,Channel.G,phi_fixed,wr0,L,sigmat2,sigmar2,Wc0,false);
end

Wc_opt = W_opt(:,1:K);
% 与SNR约束口径一致：雷达采用全部雷达流功率
wr_opt = W_opt(:,K+1:end);
end
