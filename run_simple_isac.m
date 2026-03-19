function run_simple_isac()
% 主函数：生成空间方向图(Beampattern)与帕累托前沿(Pareto Frontier)的四主方案对比
% 涵盖方案：1.RIS-Joint 2.RIS-Fixed 3.RIS-Random 4.No-RIS

fpath = which(mfilename);
if isempty(fpath), baseDir = pwd; else baseDir = fileparts(fpath); end
utilsDir = fullfile(baseDir, 'utils');
if isempty(which('evaluate_sr_baselines')), addpath(utilsDir); end

outRoot = fullfile(baseDir,'out');
if ~exist(outRoot,'dir'), mkdir(outRoot); end
runStamp = datestr(now,'yyyymmdd_HHMMSS');
outDir = fullfile(outRoot,['run_simple_isac_' runStamp]);
if ~exist(outDir,'dir'), mkdir(outDir); end

rng(1);

% 公共系统参数
M = 6; K = 4; N = 64; L = 1024;
sigmar2 = 1e-12; sigmak2 = 1e-12; sigmat2 = 1;
P_dBm = 30; P = 10.^((P_dBm-30)/10);
gammat_fixed = 10^(7/10); % 7dB 线性值

baseline = generate_baseline(M,K);
ris = generate_ris_parts(N,baseline);
Channel.hdt = baseline.hdt; Channel.Hu = baseline.Hu;
Channel.hrt = ris.hrt; Channel.G = ris.G; Channel.Hru = ris.Hru;

%% 1. 雷达空间方向图 (Radar Spatial Beampattern)
fprintf('Computing Spatial Beampattern profiles...\n');

% 四个方案相位配置与赋形求解
phi_fixed = compute_phi('fixed_target',Channel);
phi_random = compute_phi('random',Channel);

fprintf('  > [1/4] 优化 RIS-Joint 波束与相位 (最耗时)...\n');
[phi_joint,~,Wr_joint] = run_joint_snr_optimization(Channel,P,K,L,sigmat2,sigmar2,sigmak2,gammat_fixed);

fprintf('  > [2/4] 优化 RIS-Fixed 波束...\n');
[Wc_fixed, wr_fixed] = optimize_w_for_fixed_phi(Channel,phi_fixed,P,K,L,sigmat2,sigmar2,sigmak2,gammat_fixed);

fprintf('  > [3/4] 优化 RIS-Random 波束...\n');
[Wc_random, wr_random] = optimize_w_for_fixed_phi(Channel,phi_random,P,K,L,sigmat2,sigmar2,sigmak2,gammat_fixed);

Channel_noris.hdt = baseline.hdt; Channel_noris.Hu = baseline.Hu;
Channel_noris.hrt = 0; Channel_noris.G = zeros(1,M); Channel_noris.Hru = zeros(K,1);

fprintf('  > [4/4] 优化 No-RIS 波束...\n');
[Wc_noris, wr_noris] = optimize_w_for_fixed_phi(Channel_noris,1,P,K,L,sigmat2,sigmar2,sigmak2,gammat_fixed);

fprintf('  > 合成四方案对应方向图映射...\n');

% Calculate beampattern
thetaScanDeg = -90:0.2:90;
thetaScanRad = thetaScanDeg*pi/180;
aScanBS = exp(-1j*(0:M-1)'*pi*sin(thetaScanRad))/sqrt(M);
aScanRIS = exp(-1j*(0:N-1)'*pi*sin(thetaScanRad))/sqrt(N);

gain_joint = zeros(size(thetaScanDeg));
gain_fixed = zeros(size(thetaScanDeg));
gain_random = zeros(size(thetaScanDeg));
gain_noris = zeros(size(thetaScanDeg));

ampBS = norm(Channel.hdt); ampRIS = norm(Channel.hrt);

for t = 1:numel(thetaScanDeg)
    h_dir_t = ampBS*aScanBS(:,t);
    h_ris_t = ampRIS*aScanRIS(:,t);
    
    h_eff_joint_t = h_dir_t.' + h_ris_t.'*diag(phi_joint)*Channel.G;
    h_eff_fixed_t = h_dir_t.' + h_ris_t.'*diag(phi_fixed)*Channel.G;
    h_eff_random_t = h_dir_t.' + h_ris_t.'*diag(phi_random)*Channel.G;

    gain_joint(t) = norm(h_eff_joint_t*Wr_joint)^2;
    gain_fixed(t) = norm(h_eff_fixed_t*[Wc_fixed, wr_fixed])^2;
    gain_random(t) = norm(h_eff_random_t*[Wc_random, wr_random])^2;
    gain_noris(t) = norm(h_dir_t.'*[Wc_noris, wr_noris])^2;
end

gainRef = max(gain_joint(:));
gainDb_joint = 10*log10(gain_joint./(gainRef+1e-12) + 1e-12);
gainDb_fixed = 10*log10(gain_fixed./(gainRef+1e-12) + 1e-12);
gainDb_random = 10*log10(gain_random./(gainRef+1e-12) + 1e-12);
gainDb_noris = 10*log10(gain_noris./(gainRef+1e-12) + 1e-12);

targetAngleDeg = estimate_doa_from_channel(Channel.hdt);
userAngleDeg = zeros(K,1);
for k=1:K
    userAngleDeg(k) = estimate_doa_from_channel(Channel.Hu(k,:).');
end

fig_bp = figure('Color','w');
plot(thetaScanDeg,gainDb_joint,'-','LineWidth',1.5,'Color',[0.85 0.15 0.15]); hold on;
plot(thetaScanDeg,gainDb_fixed,'--','LineWidth',1.5,'Color',[0.90 0.50 0.10]);
plot(thetaScanDeg,gainDb_random,'-.','LineWidth',1.5,'Color',[0.10 0.55 0.10]);
plot(thetaScanDeg,gainDb_noris,':','LineWidth',1.5,'Color',[0.20 0.20 0.20]);
grid on;
xlabel('Angle (deg)'); ylabel('Equivalent Beampattern Gain (dB)');
xlim([-90 90]); ylim([-60 5]);
ylbp = ylim;
plot([targetAngleDeg targetAngleDeg], ylbp, ':', 'LineWidth',1.2,'Color',[0.2 0.2 0.2]);
for k = 1:K
    plot([userAngleDeg(k) userAngleDeg(k)], ylbp, '--', 'LineWidth',0.8,'Color',[0.6 0.6 0.6]);
end
text(targetAngleDeg, ylbp(2)-3, 'Target', 'HorizontalAlignment','center');
legend('RIS-Joint','RIS-Fixed','RIS-Random','No-RIS','Location','southoutside');
title('Radar Spatial Beampattern (Optimized)');
print(fig_bp, fullfile(outDir,'radar_spatial_beampattern.png'), '-dpng','-r300');
savefig(fig_bp, fullfile(outDir,'radar_spatial_beampattern.fig'));


%% 2. 帕累托前沿图 (Pareto Frontier)
fprintf('Computing Pareto Frontier limits...\n');
gammat_pf_dB = 0:1:15;
SR_pf_joint = zeros(size(gammat_pf_dB)); SNR_pf_joint = zeros(size(gammat_pf_dB));
SR_pf_fixed = zeros(size(gammat_pf_dB)); SNR_pf_fixed = zeros(size(gammat_pf_dB));
SR_pf_random = zeros(size(gammat_pf_dB)); SNR_pf_random = zeros(size(gammat_pf_dB));
SR_pf_noris = zeros(size(gammat_pf_dB)); SNR_pf_noris = zeros(size(gammat_pf_dB));

prev_state = [];
for i = 1:numel(gammat_pf_dB)
    gt = 10^(gammat_pf_dB(i)/10);
    fprintf('[%d/%d] 正在计算 Pareto 边界点 Gamma_t = %g dB...\n', i, numel(gammat_pf_dB), gammat_pf_dB(i));
    
    [out_i, prev_state] = evaluate_sr_baselines(Channel,baseline,P,K,L,sigmat2,sigmar2,sigmak2,gt,prev_state);
    
    SR_pf_joint(i) = out_i.rate.ris_joint;
    SNR_pf_joint(i) = out_i.snr.ris_joint;
    
    SR_pf_fixed(i) = out_i.rate.ris_fixed;
    SNR_pf_fixed(i) = out_i.snr.ris_fixed;
    
    SR_pf_random(i) = out_i.rate.ris_random;
    SNR_pf_random(i) = out_i.snr.ris_random;
    
    SR_pf_noris(i) = out_i.rate.no_ris;
    SNR_pf_noris(i) = out_i.snr.no_ris;
    
    fprintf('  -> Joint SR: %.4f | Fixed SR: %.4f | Random SR: %.4f | No-RIS SR: %.4f\n', ...
            out_i.rate.ris_joint, out_i.rate.ris_fixed, out_i.rate.ris_random, out_i.rate.no_ris);
end

fig_pf = figure('Color','w');
plot(10*log10(SNR_pf_joint),SR_pf_joint,'-o','LineWidth',2.0,'Color',[0.85 0.15 0.15]); hold on;
plot(10*log10(SNR_pf_fixed),SR_pf_fixed,'-s','LineWidth',2.0,'Color',[0.90 0.50 0.10]);
plot(10*log10(SNR_pf_random),SR_pf_random,'-^','LineWidth',2.0,'Color',[0.10 0.55 0.10]);
plot(10*log10(SNR_pf_noris),SR_pf_noris,'-d','LineWidth',2.0,'Color',[0.20 0.20 0.20]);
grid on; hold off;
xlabel('Radar SNR (dB)'); ylabel('Sum-rate (bps/Hz)');
legend('RIS-Joint frontier','RIS-Fixed frontier','RIS-Random frontier','No-RIS frontier','Location','best');
title('ISAC Pareto Frontier (Optimal, sweep \Gamma_r)');
print(fig_pf, fullfile(outDir,'pareto_frontier.png'), '-dpng','-r300');
savefig(fig_pf, fullfile(outDir,'pareto_frontier.fig'));

fprintf('Run completed successfully. Results saved in %s\n', outDir);
end

function angDeg = estimate_doa_from_channel(h)
h = h(:);
if numel(h) < 2
    angDeg = 0; return;
end
phaseDiff = angle(h(2)/h(1));
s = -phaseDiff/pi;
s = max(min(real(s),1),-1);
angDeg = asind(s);
end