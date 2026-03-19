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
    % 动态计算稳健的初始 eta0：根据目标门限估算雷达功率，并给一点余量确保初值可行
    h_eff_t_norm2 = norm(Channel.hdt.' + Channel.hrt.'*diag(phi_fixed)*Channel.G)^2;
    % 根据 SNR 公式推导需要的 P_rad
    Prad_req = (gammat * sigmar2) / (max(h_eff_t_norm2, 1e-15) * L * sigmat2);
    eta0 = (Prad_req / P) * 1.5; % 给 1.5 倍的功率余量，确认为可行解开始优化
    eta0 = max(0.05, min(0.9, eta0)); % 底层保持 5% 雷达功率，防止梯度消失导致 NaN
    
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
    fprintf('\n=> [FATAL ERROR] optimize_w_for_fixed_phi 抛出异常: %s\n', ME.message);
    fprintf('=>        [Stack] File: %s, Line: %d\n\n', ME.stack(1).file, ME.stack(1).line);
    W_opt = W0;
    Vsr_opt = sum_rate(Hk0,Wc0,sigmak2,wr0,false);
    gammat_opt = radar_snr(Channel.hdt,Channel.hrt,Channel.G,phi_fixed,wr0,L,sigmat2,sigmar2,Wc0,false);
end

Wc_opt = W_opt(:,1:K);
% 与SNR约束口径一致：雷达采用全部雷达流功率
wr_opt = W_opt(:,K+1:end);
end
