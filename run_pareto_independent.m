function run_pareto_independent()
% run_pareto_independent - 独立点评估帕累托前沿曲线
% 
% 核心逻辑：
% 1. 雷达要求门限 Gamma_t 每隔 2 dB 计算一次。
% 2. 取消了门限点之间的温启动（即 6dB 不会拿 4dB 的结果去当初始值，杜绝跨点干扰）。
% 3. 在任意单个考点内：
%      - RIS-Fixed、RIS-Random、No-RIS: 每次均使用最原生的均匀功率分配法 (eta=0.5) 独立开荒寻找波束初值。
%      - RIS-Joint (联合优化): 每次严格白嫖旁边刚刚算出的 RIS-Fixed 最优相位和波束作为完美垫脚石，大幅缩减自己单独找初始可行解的时间。

fpath = which(mfilename);
if isempty(fpath), baseDir = pwd; else baseDir = fileparts(fpath); end
utilsDir = fullfile(baseDir, 'utils');
refDir = fullfile(baseDir, 'reference');
try rmpath(refDir); catch; end
clear functions; % 强制刷新缓存，识别已删除的幽灵文件
if isempty(which('evaluate_sr_baselines')), addpath(utilsDir); end

outRoot = fullfile(baseDir,'out');
if ~exist(outRoot,'dir'), mkdir(outRoot); end
runStamp = datestr(now,'yyyymmdd_HHMMSS');
outDir = fullfile(outRoot,['run_pareto_independent_' runStamp]);
if ~exist(outDir,'dir'), mkdir(outDir); end

rng(1);

% 公共系统参数
M = 6; K = 4; N = 64; L = 1024;
sigmar2 = 1e-12; sigmak2 = 1e-12; sigmat2 = 1;
P_dBm = 30; P = 10.^((P_dBm-30)/10);

baseline = generate_baseline(M,K);
ris = generate_ris_parts(N,baseline);
Channel.hdt = baseline.hdt; Channel.Hu = baseline.Hu;
Channel.hrt = ris.hrt; Channel.G = ris.G; Channel.Hru = ris.Hru;

%% ===== 第一阶段：评估纯通信天花板 (Comm-only Ceiling) =====
% 作为后续低要求门限的“完美起跑器”
fprintf('\n=======================================================\n');
fprintf('>>> 1. 正在榨取绝对纯通信性能上限 (Comm-only Baseline)... \n');
Prms.M = M; Prms.N = N; Prms.K = K; Prms.sigmar2 = sigmar2; 
Prms.sigmak2 = sigmak2; Prms.sigmat2 = sigmat2; Prms.Nmax = 100; 
Prms.res_th = 5e-4; Prms.L = L; Prms.P = P;

phi0_comm = compute_phi('random',Channel);
Hk0_comm = Channel.Hu + Channel.Hru*diag(phi0_comm)*Channel.G;
W0_comm = Hk0_comm'/(Hk0_comm*Hk0_comm');
for k = 1:K, W0_comm(:,k) = W0_comm(:,k)/norm(W0_comm(:,k)); end
W0_comm = W0_comm * sqrt(P/K);

[W_comm,phi_comm,Vsr_comm,gammat_comm] = get_SR_comm_only(Prms,Channel,phi0_comm,W0_comm);
SR_comm = Vsr_comm(end);
SNR_comm_dB = 10*log10(abs(gammat_comm)+1e-12);
fprintf('  -> 🌟 纯通信极限 (Comm-only): SR: %.4f, 自然雷达 SNR: %.2f dB\n\n', SR_comm, SNR_comm_dB);

% 预装配完美起跑器（专供低于自然 SNR 的考核点使用）
prev_state_comm = struct();
prev_state_comm.phi_joint = phi_comm;
% W_comm 只有通信维度的 M x K，需扩充全零的雷达波束维度 M x M 以兼容 Joint 优化器的入参尺寸
prev_state_comm.W_joint = [W_comm, zeros(M, M)];


%% ===== 第二阶段：计算帕累托前沿 (Pareto Frontier) =====
fprintf('\n=======================================================\n');
fprintf('>>> 2. 开始独立点 Pareto 扫描 (扫描范围 0:2:14 dB)\n');
fprintf('=======================================================\n\n');

gammat_pf_dB = 0:2:14; % 定制需求：每隔 2 dB 画一次

SR_pf_joint = zeros(size(gammat_pf_dB)); SNR_pf_joint = zeros(size(gammat_pf_dB));
SR_pf_fixed = zeros(size(gammat_pf_dB)); SNR_pf_fixed = zeros(size(gammat_pf_dB));
SR_pf_random = zeros(size(gammat_pf_dB)); SNR_pf_random = zeros(size(gammat_pf_dB));
SR_pf_noris = zeros(size(gammat_pf_dB)); SNR_pf_noris = zeros(size(gammat_pf_dB));

% Joint 的温启动链：第一次用 Comm-only，后续每次继承上一轮 Joint 的输出
prev_state_joint = prev_state_comm;  % 种子 = 纯通信天花板

for i = 1:numel(gammat_pf_dB)
    gt = 10^(gammat_pf_dB(i)/10);
    fprintf('>>> [%d/%d] 正在评估独立目标考点 Gamma_t = %g dB...\n', i, numel(gammat_pf_dB), gammat_pf_dB(i));
    
    % Joint 方案：始终传入上一轮的结果做温启动（第一轮 = Comm-only）
    fprintf('  * Joint 初始化来源: 上一轮 Joint 输出 (首轮=Comm-only)\n');
    [out_i, cur_state] = evaluate_sr_baselines(Channel,baseline,P,K,L,sigmat2,sigmar2,sigmak2,gt,prev_state_joint);
    
    % 更新温启动链：将本轮 Joint 的输出传递给下一轮
    if out_i.feasible.ris_joint
        prev_state_joint = struct();
        prev_state_joint.phi_joint = cur_state.phi_joint;
        prev_state_joint.W_joint = cur_state.W_joint;
    end
    % 如果本轮 Joint 不可行，保持上一轮的 prev_state 不变（自动降级）
    
    SR_pf_joint(i) = out_i.rate.ris_joint; SNR_pf_joint(i) = out_i.snr.ris_joint;
    SR_pf_fixed(i) = out_i.rate.ris_fixed; SNR_pf_fixed(i) = out_i.snr.ris_fixed;
    SR_pf_random(i) = out_i.rate.ris_random; SNR_pf_random(i) = out_i.snr.ris_random;
    SR_pf_noris(i) = out_i.rate.no_ris; SNR_pf_noris(i) = out_i.snr.no_ris;
    
    fprintf('  -> 本轮结果: Joint SR: %.4f | Fixed SR: %.4f | Random SR: %.4f | No-RIS SR: %.4f\n\n', ...
            out_i.rate.ris_joint, out_i.rate.ris_fixed, out_i.rate.ris_random, out_i.rate.no_ris);
end

% 绘图逻辑
fig_pf = figure('Color','w');
plot(10*log10(SNR_pf_joint),SR_pf_joint,'-o','LineWidth',2.0,'Color',[0.85 0.15 0.15]); hold on;
plot(10*log10(SNR_pf_fixed),SR_pf_fixed,'-s','LineWidth',2.0,'Color',[0.90 0.50 0.10]);
plot(10*log10(SNR_pf_random),SR_pf_random,'-^','LineWidth',2.0,'Color',[0.10 0.55 0.10]);
plot(10*log10(SNR_pf_noris),SR_pf_noris,'-d','LineWidth',2.0,'Color',[0.20 0.20 0.20]);
% 标记纯通信物理极限
plot(SNR_comm_dB, SR_comm, 'p', 'MarkerSize', 18, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k');
xline(SNR_comm_dB, '--k', 'Comm-only Natural SNR', 'LabelVerticalAlignment', 'bottom');
yline(SR_comm, '--k', 'Comm-only SR Ceiling', 'LabelHorizontalAlignment', 'left');
grid on; hold off;

xlabel('Radar SNR \Gamma_t (dB)'); ylabel('Sum-rate (bps/Hz)');
legend('RIS-Joint frontier','RIS-Fixed frontier','RIS-Random frontier','No-RIS frontier','Location','best');
title('ISAC Pareto Frontier (Independent Evaluation)');

print(fig_pf, fullfile(outDir,'pareto_frontier_independent.png'), '-dpng','-r300');
savefig(fig_pf, fullfile(outDir,'pareto_frontier_independent.fig'));

fprintf('运行圆满完成！独立图表及数据已存入: %s\n', outDir);

end
