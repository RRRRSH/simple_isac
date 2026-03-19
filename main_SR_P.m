function main_SR_P(gammatFixed_dB,outDir)
% main_SR_P - 绘制 Sum-rate 随发射功率 P_t 变化曲线（四主基线对比）

fpath = which(mfilename);
if isempty(fpath), baseDir = pwd; else baseDir = fileparts(fpath); end
utilsDir = fullfile(baseDir,'utils');
if isempty(which('evaluate_sr_baselines')), addpath(utilsDir); end

rng(1);

if nargin < 1 || isempty(gammatFixed_dB)
    gammatFixed_dB = 7;
end
if nargin < 2 || isempty(outDir)
    outRoot = fullfile(baseDir,'out');
    if ~exist(outRoot,'dir'), mkdir(outRoot); end
    runStamp = datestr(now,'yyyymmdd_HHMMSS');
    outDir = fullfile(outRoot,['main_SR_P_' runStamp]);
end
if ~exist(outDir,'dir'), mkdir(outDir); end
fprintf('Results will be saved to: %s\n', outDir);

M = 6; K = 4; N = 64; L = 1024;
sigmar2 = 1e-12; sigmak2 = 1e-12; sigmat2 = 1;

Pt_dBm = 20:2:30;
Pt_W = 10.^((Pt_dBm-30)/10);
gammat = 10^(gammatFixed_dB/10);

baseline = generate_baseline(M,K);
ris = generate_ris_parts(N,baseline);
Channel.hdt = baseline.hdt; Channel.Hu = baseline.Hu;
Channel.hrt = ris.hrt; Channel.G = ris.G; Channel.Hru = ris.Hru;

SR_no_ris = zeros(size(Pt_W));
SR_ris_random = zeros(size(Pt_W));
SR_ris_fixed = zeros(size(Pt_W));
SR_ris_joint = zeros(size(Pt_W));

prev_state = [];
for i = 1:numel(Pt_W)
    [out_i, prev_state] = evaluate_sr_baselines(Channel,baseline,Pt_W(i),K,L,sigmat2,sigmar2,sigmak2,gammat,prev_state);
    SR_no_ris(i) = out_i.rate.no_ris;
    SR_ris_random(i) = out_i.rate.ris_random;
    SR_ris_fixed(i) = out_i.rate.ris_fixed;
    SR_ris_joint(i) = out_i.rate.ris_joint;
    
    fprintf('[%d/%d] Pt=%.1f dBm, No-RIS SR=%.4f, Random SR=%.4f, Fixed SR=%.4f, Joint SR=%.4f\n', ...
        i, numel(Pt_W), Pt_dBm(i), SR_no_ris(i), SR_ris_random(i), SR_ris_fixed(i), SR_ris_joint(i));
end

fig = figure('Color','w');
plot(Pt_dBm,SR_ris_joint,'-o','LineWidth',1.8,'Color',[0.85 0.15 0.15]); hold on;
plot(Pt_dBm,SR_ris_fixed,'-s','LineWidth',1.8,'Color',[0.90 0.50 0.10]);
plot(Pt_dBm,SR_ris_random,'-^','LineWidth',1.8,'Color',[0.10 0.55 0.10]);
plot(Pt_dBm,SR_no_ris,'-d','LineWidth',1.8,'Color',[0.20 0.20 0.20]);
grid on; hold off;

xlabel('Transmit power P_t (dBm)');
ylabel('Sum-rate (bps/Hz)');
legend('RIS-Joint (Proposed)','RIS-Fixed','RIS-Random','No-RIS','Location','best');
title(sprintf('Sum-rate vs Transmit Power (\\Gamma_t = %g dB)', gammatFixed_dB));

print(fig, fullfile(outDir,'sumrate_vs_Pt.png'), '-dpng','-r300');
savefig(fig, fullfile(outDir,'sumrate_vs_Pt.fig'));

save(fullfile(outDir,'sumrate_vs_Pt_data.mat'), ...
    'Pt_dBm','Pt_W','gammatFixed_dB','gammat', ...
    'SR_no_ris','SR_ris_random','SR_ris_fixed','SR_ris_joint');
end
