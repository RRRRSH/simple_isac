function main_SR_N(gammatFixed_dB,runMode,outDir)
% main_SR_N - 绘制 Sum-rate 随 RIS 单元数 N 变化曲线（四方案对比）
% 可选输入:
%   gammatFixed_dB: 固定门限图使用的雷达门限(dB)，默认 7
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

if nargin < 1 || isempty(gammatFixed_dB)
    gammatFixed_dB = 7;
end
if nargin < 2 || isempty(runMode)
    runMode = 'all';
end

if nargin < 3 || isempty(outDir)
    outRoot = fullfile(baseDir,'out');
    if ~exist(outRoot,'dir'), mkdir(outRoot); end
    runStamp = datestr(now,'yyyymmdd_HHMMSS');
    outDir = fullfile(outRoot,['main_SR_N_' runStamp]);
end
if ~exist(outDir,'dir'), mkdir(outDir); end
fprintf('Results will be saved to: %s\n', outDir);
fprintf('Run mode: %s, fixed threshold: %.3f dB\n', char(runMode), gammatFixed_dB);
isNoRISOnly = any(lower(string(runMode)) == ["noris_only","no_ris_only"]);
isNoRISVsFixed = any(lower(string(runMode)) == ["noris_vs_fixed","no_ris_vs_fixed"]);
isNoRISVsFixedNSPMRT = any(lower(string(runMode)) == ["noris_vs_fixed_nsp_vs_mrt","no_ris_vs_fixed_nsp_vs_mrt"]);

M = 6;
K = 4;
L = 1024;
sigmar2 = 1e-12;
sigmak2 = 1e-12;
sigmat2 = 1;
includeRadarInterferenceInRate = false;
includeCommInterferenceInRadar = false;

N_list = [16 25 36 49 64];
P_dBm = 30;
P = 10.^((P_dBm-30)/10);

baseline = generate_baseline(M,K);
Nmax = max(N_list);
ris_master = generate_ris_parts(Nmax,baseline);

% 固定门限图使用固定阈值（默认7 dB，可通过入参修改）
gammat = 10^(gammatFixed_dB/10);

SR_proposed = zeros(size(N_list));
SR_separate = zeros(size(N_list));
SR_comm_only = zeros(size(N_list));
SR_no_ris = zeros(size(N_list));
SR_fixed_nsp = zeros(size(N_list));
SR_fixed_mrt = zeros(size(N_list));

for i = 1:numel(N_list)
    N = N_list(i);
    ris = generate_ris_parts(N,baseline,ris_master);
    Channel.hdt = baseline.hdt;
    Channel.Hu = baseline.Hu;
    Channel.hrt = ris.hrt;
    Channel.G = ris.G;
    Channel.Hru = ris.Hru;

    out_i = evaluate_sr_baselines(Channel,baseline,P,K,L,sigmat2,sigmar2,sigmak2,gammat, ...
        includeRadarInterferenceInRate,includeCommInterferenceInRadar,runMode);
    SR_proposed(i) = out_i.rate.proposed;
    SR_separate(i) = out_i.rate.separate;
    SR_comm_only(i) = out_i.rate.comm_only;
    SR_no_ris(i) = out_i.rate.no_ris;
    SR_fixed_nsp(i) = out_i.rate.fixed_nsp;
    SR_fixed_mrt(i) = out_i.rate.fixed_mrt;
    if isNoRISVsFixedNSPMRT
        fprintf('[%d/%d] N=%d, No-RIS SR=%g, RIS+MRT SR=%g, RIS+NSP SR=%g\n', ...
            i, numel(N_list), N, SR_no_ris(i), SR_fixed_mrt(i), SR_fixed_nsp(i));
    elseif isNoRISVsFixed
        fprintf('[%d/%d] N=%d, No-RIS SR=%g, RIS-Fixed SR=%g\n', ...
            i, numel(N_list), N, SR_no_ris(i), SR_separate(i));
    else
        fprintf('[%d/%d] N=%d, No-RIS SR=%g\n', i, numel(N_list), N, SR_no_ris(i));
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
    plot(N_list,SR_no_ris_plot,'-d','LineWidth',1.9,'Color',[0.20 0.20 0.20]); hold on;
elseif isNoRISVsFixed
    plot(N_list,SR_fixed_plot,'-s','LineWidth',1.8,'Color',[0.85 0.33 0.10]); hold on;
    plot(N_list,SR_no_ris_plot,'-d','LineWidth',1.8,'Color',[0.20 0.20 0.20]);
elseif isNoRISVsFixedNSPMRT
    plot(N_list,SR_fixed_mrt_plot,'-^','LineWidth',1.8,'Color',[0.95 0.50 0.10]); hold on;
    plot(N_list,SR_fixed_nsp_plot,'-s','LineWidth',1.8,'Color',[0.10 0.55 0.10]);
    plot(N_list,SR_no_ris_plot,'-d','LineWidth',1.8,'Color',[0.20 0.20 0.20]);
else
    plot(N_list,SR_proposed,'-o','LineWidth',1.7,'Color',[0.00 0.45 0.74]); hold on;
    plot(N_list,SR_separate,'-s','LineWidth',1.7,'Color',[0.85 0.33 0.10]);
    plot(N_list,SR_comm_only,'-^','LineWidth',1.7,'Color',[0.47 0.67 0.19]);
    plot(N_list,SR_no_ris_plot,'-d','LineWidth',1.7,'Color',[0.20 0.20 0.20]);
end
hold off;
grid on;
xlabel('Number of RIS elements N');
ylabel('Sum-rate (bps/Hz)');
if isNoRISOnly
    legend('No-RIS / BF-only','Location','best');
    title(sprintf('Sum-rate vs RIS Elements (No-RIS only, P_t = %d dBm)', P_dBm));
elseif isNoRISVsFixed
    legend('RIS-Fixed','No-RIS / BF-only','Location','best');
    title(sprintf('Sum-rate vs RIS Elements (No-RIS vs RIS-Fixed, P_t = %d dBm)', P_dBm));
elseif isNoRISVsFixedNSPMRT
    legend('RIS+MRT','RIS+NSP','No-RIS / BF-only','Location','best');
    title(sprintf('Sum-rate vs RIS Elements (No-RIS vs RIS+MRT vs RIS+NSP, P_t = %d dBm)', P_dBm));
else
    legend('Proposed (Joint)','Separate','Comm-only','No-RIS / BF-only','Location','best');
    title(sprintf('Sum-rate vs RIS Elements (P_t = %d dBm)', P_dBm));
end

print(fig, fullfile(outDir,'sumrate_vs_N.png'), '-dpng','-r300');
savefig(fig, fullfile(outDir,'sumrate_vs_N.fig'));

save(fullfile(outDir,'sumrate_vs_N_data.mat'), ...
    'N_list','P_dBm','P','gammatFixed_dB','gammat','runMode', ...
    'SR_proposed','SR_separate','SR_comm_only','SR_no_ris','SR_fixed_nsp','SR_fixed_mrt', ...
    'SR_no_ris_plot','SR_fixed_plot','SR_fixed_nsp_plot','SR_fixed_mrt_plot', ...
    'infeasibleMaskNoRIS','infeasibleMaskFixed','infeasibleMaskFixedNSP','infeasibleMaskFixedMRT');
end
