function run_noris_only()
% run_noris_only
% 快速仅运行 No-RIS 场景：
% 1) 固定基站总功率 P=1W
% 2) 在满足雷达SNR门限下优化 No-RIS 波束 (P1在No-RIS下的简化)
% 3) 在终端输出清晰进度日志

baseDir = fileparts(mfilename('fullpath'));
if isempty(baseDir), baseDir = pwd; end
addpath(fullfile(baseDir,'utils'));

fprintf('\n========== No-RIS Quick Run ==========%s\n', newline);
fprintf('[1/6] 初始化参数...\n');
rng(1);
M = 6;
K = 4;
L = 1024;
sigmar2 = 1e-12;
sigmak2 = 1e-12;
sigmat2 = 1;
P = 1.0;  % 固定1W

fprintf('      M=%d, K=%d, L=%d, P=%.3f W\n', M, K, L, P);

fprintf('[2/6] 生成基线信道...\n');
baseline = generate_baseline(M,K);

fprintf('[3/6] 构造No-RIS信道...\n');
% 通过退化信道实现No-RIS：Hru*diag(phi)*G = 0
Channel.hdt = baseline.hdt;
Channel.Hu = baseline.Hu;
Channel.hrt = 0;            % length=1，兼容优化器
Channel.G = zeros(1,M);     % 1xM
Channel.Hru = zeros(K,1);   % Kx1
phi_noris = 1;

fprintf('[4/6] 设置雷达SNR门限...\n');
% 统一雷达SNR门限：固定为7 dB（线性值）
gammat = 10^(7/10);
fprintf('      gammat(linear)=%.3e, gammat(dB)=%.3f dB\n', gammat, 10*log10(gammat));

fprintf('[5/6] 运行No-RIS波束优化（终端会输出迭代进度）...\n');
[Wc_opt,wr_opt,~,Vsr_opt,gammat_opt] = optimize_w_for_fixed_phi( ...
    Channel,phi_noris,P,K,L,sigmat2,sigmar2,sigmak2,gammat,true);

fprintf('[6/6] 计算并汇总性能...\n');
Hk_noris = Channel.Hu;
SR = sum_rate(Hk_noris,Wc_opt,sigmak2,wr_opt,false);
SNR_legacy = radar_snr_noris(Channel.hdt,wr_opt,L,sigmat2,sigmar2,Wc_opt,false);
ptotal = norm(Wc_opt,'fro')^2 + norm(wr_opt,2)^2;

fprintf('\n===== No-RIS Result Summary =====\n');
fprintf('Sum-rate              : %.6f bps/Hz\n', SR);
fprintf('Radar SNR (optimizer) : %.6f\n', gammat_opt);
fprintf('SNR threshold met?    : %d\n', gammat_opt >= gammat);
fprintf('Radar SNR (legacy fn) : %.3e (reference only)\n', SNR_legacy);
fprintf('Total TX power used   : %.6f W\n', ptotal);
fprintf('Power constraint met? : %d\n', ptotal <= P + 1e-6);
if ~isempty(Vsr_opt)
    fprintf('Final optimizer Vsr   : %.6f\n', Vsr_opt(end));
end
fprintf('Final optimizer SNR   : %.6f\n', gammat_opt);
fprintf('=================================\n\n');
end
