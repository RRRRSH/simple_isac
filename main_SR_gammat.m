function main_SR_gammat(gammat_dB,runMode,outDir)
% main_SR_gammat - 绘制 Sum-rate 随雷达门限 Gamma_t 变化曲线（四方案对比）
% 可选输入:
%   gammat_dB: 门限扫描向量(dB)，默认 0:2:10
%   runMode: 'all'(默认) / 'noris_only' / 'noris_vs_fixed' / 'noris_vs_fixed_nsp_vs_mrt'
%   outDir: 输出目录。未提供时自动按时间戳创建

fpath = which(mfilename);
if isempty(fpath)
    baseDir = pwd;
else
    baseDir = fileparts(fpath);
end
utilsDir = fullfile(baseDir,'utils');
if isempty(which('evaluate_sr_baselines'))
    addpath(utilsDir);
end

rng(1);

if nargin < 1 || isempty(gammat_dB)
    gammat_dB = 0:2:10;
end
if nargin < 2 || isempty(runMode)
    runMode = 'all';
end

if nargin < 3 || isempty(outDir)
    outRoot = fullfile(baseDir,'out');
    if ~exist(outRoot,'dir'), mkdir(outRoot); end
    runStamp = datestr(now,'yyyymmdd_HHMMSS');
    outDir = fullfile(outRoot,['main_SR_gammat_' runStamp]);
end
if ~exist(outDir,'dir'), mkdir(outDir); end
fprintf('Results will be saved to: %s\n', outDir);
fprintf('Run mode: %s, threshold sweep(dB): [%s]\n', char(runMode), num2str(gammat_dB));
isNoRISOnly = any(lower(string(runMode)) == ["noris_only","no_ris_only"]);
isNoRISVsFixed = any(lower(string(runMode)) == ["noris_vs_fixed","no_ris_vs_fixed"]);
isNoRISVsFixedNSPMRT = any(lower(string(runMode)) == ["noris_vs_fixed_nsp_vs_mrt","no_ris_vs_fixed_nsp_vs_mrt"]);

M = 6;
K = 4;
N = 64;
L = 1024;
sigmar2 = 1e-12;
sigmak2 = 1e-12;
sigmat2 = 1;
includeRadarInterferenceInRate = false;
includeCommInterferenceInRadar = false;

P_dBm = 30;
P = 10.^((P_dBm-30)/10);

baseline = generate_baseline(M,K);
ris = generate_ris_parts(N,baseline);
Channel.hdt = baseline.hdt;
Channel.Hu = baseline.Hu;
Channel.hrt = ris.hrt;
Channel.G = ris.G;
Channel.Hru = ris.Hru;

% 统一门限口径：不做缩放，直接使用线性值
gammat_linear = 10.^(gammat_dB/10);

SR_proposed = zeros(size(gammat_dB));
SR_separate = zeros(size(gammat_dB));
SR_comm_only = zeros(size(gammat_dB));
SR_no_ris = zeros(size(gammat_dB));
SR_fixed_nsp = zeros(size(gammat_dB));
SR_fixed_mrt = zeros(size(gammat_dB));

for i = 1:numel(gammat_dB)
    out_i = evaluate_sr_baselines(Channel,baseline,P,K,L,sigmat2,sigmar2,sigmak2,gammat_linear(i), ...
        includeRadarInterferenceInRate,includeCommInterferenceInRadar,runMode);
    SR_proposed(i) = out_i.rate.proposed;
    SR_separate(i) = out_i.rate.separate;
    SR_comm_only(i) = out_i.rate.comm_only;
    SR_no_ris(i) = out_i.rate.no_ris;
    SR_fixed_nsp(i) = out_i.rate.fixed_nsp;
    SR_fixed_mrt(i) = out_i.rate.fixed_mrt;
    if isNoRISVsFixedNSPMRT
        fprintf('[%d/%d] gamma=%.2f dB, No-RIS SR=%g, RIS+MRT SR=%g, RIS+NSP SR=%g\n', ...
            i, numel(gammat_dB), gammat_dB(i), SR_no_ris(i), SR_fixed_mrt(i), SR_fixed_nsp(i));
    elseif isNoRISVsFixed
        fprintf('[%d/%d] gamma=%.2f dB, No-RIS SR=%g, RIS-Fixed SR=%g\n', ...
            i, numel(gammat_dB), gammat_dB(i), SR_no_ris(i), SR_separate(i));
    else
        fprintf('[%d/%d] gamma=%.2f dB, No-RIS SR=%g\n', i, numel(gammat_dB), gammat_dB(i), SR_no_ris(i));
    end
end

infeasibleMaskNoRIS = isnan(SR_no_ris);
SR_no_ris_plot = SR_no_ris;
SR_no_ris_plot(infeasibleMaskNoRIS) = 0;
infeasibleMaskFixed = isnan(SR_separate);
SR_fixed_plot = SR_separate;
SR_fixed_plot(infeasibleMaskFixed) = 0;
infeasibleMaskFixedNSP = isnan(SR_fixed_nsp);
SR_fixed_nsp_plot = SR_fixed_nsp;
SR_fixed_nsp_plot(infeasibleMaskFixedNSP) = 0;
infeasibleMaskFixedMRT = isnan(SR_fixed_mrt);
SR_fixed_mrt_plot = SR_fixed_mrt;
SR_fixed_mrt_plot(infeasibleMaskFixedMRT) = 0;

fig = figure('Color','w');
if isNoRISOnly
    plot(gammat_dB,SR_no_ris_plot,'-d','LineWidth',1.9,'Color',[0.20 0.20 0.20]); hold on;
elseif isNoRISVsFixed
    plot(gammat_dB,SR_fixed_plot,'-s','LineWidth',1.8,'Color',[0.85 0.33 0.10]); hold on;
    plot(gammat_dB,SR_no_ris_plot,'-d','LineWidth',1.8,'Color',[0.20 0.20 0.20]);
elseif isNoRISVsFixedNSPMRT
    plot(gammat_dB,SR_fixed_mrt_plot,'-^','LineWidth',1.8,'Color',[0.95 0.50 0.10]); hold on;
    plot(gammat_dB,SR_fixed_nsp_plot,'-s','LineWidth',1.8,'Color',[0.10 0.55 0.10]);
    plot(gammat_dB,SR_no_ris_plot,'-d','LineWidth',1.8,'Color',[0.20 0.20 0.20]);
else
    plot(gammat_dB,SR_proposed,'-o','LineWidth',1.7,'Color',[0.00 0.45 0.74]); hold on;
    plot(gammat_dB,SR_separate,'-s','LineWidth',1.7,'Color',[0.85 0.33 0.10]);
    plot(gammat_dB,SR_comm_only,'-^','LineWidth',1.7,'Color',[0.47 0.67 0.19]);
    plot(gammat_dB,SR_no_ris_plot,'-d','LineWidth',1.7,'Color',[0.20 0.20 0.20]);
end
hold off;
grid on;
xlabel('Radar threshold \Gamma_t (dB)');
ylabel('Sum-rate (bps/Hz)');
if isNoRISOnly
    legend('No-RIS / BF-only','Location','best');
    title(sprintf('Communication-Sensing Trade-off (No-RIS only, P_t = %d dBm, N = %d)', P_dBm, N));
elseif isNoRISVsFixed
    legend('RIS-Fixed','No-RIS / BF-only','Location','best');
    title(sprintf('Communication-Sensing Trade-off (No-RIS vs RIS-Fixed, P_t = %d dBm, N = %d)', P_dBm, N));
elseif isNoRISVsFixedNSPMRT
    legend('RIS+MRT','RIS+NSP','No-RIS / BF-only','Location','best');
    title(sprintf('Communication-Sensing Trade-off (No-RIS vs RIS+MRT vs RIS+NSP, P_t = %d dBm, N = %d)', P_dBm, N));
else
    legend('Proposed (Joint)','Separate','Comm-only','No-RIS / BF-only','Location','best');
    title(sprintf('Communication-Sensing Trade-off (P_t = %d dBm, N = %d)', P_dBm, N));
end

print(fig, fullfile(outDir,'sumrate_vs_gammat.png'), '-dpng','-r300');
savefig(fig, fullfile(outDir,'sumrate_vs_gammat.fig'));

save(fullfile(outDir,'sumrate_vs_gammat_data.mat'), ...
    'gammat_dB','gammat_linear','P_dBm','P','N','runMode', ...
    'SR_proposed','SR_separate','SR_comm_only','SR_no_ris','SR_fixed_nsp','SR_fixed_mrt', ...
    'SR_no_ris_plot','SR_fixed_plot','SR_fixed_nsp_plot','SR_fixed_mrt_plot', ...
    'infeasibleMaskNoRIS','infeasibleMaskFixed','infeasibleMaskFixedNSP','infeasibleMaskFixedMRT');
end
