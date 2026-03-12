function run_simple_isac(varargin)
% 主函数：运行RIS辅助的ISAC(集成感知与通信)系统仿真
% 功能：对比不同RIS相位配置策略下的通信速率和雷达性能
% 用法：
%   run_simple_isac()                 % 默认同时启用两种模式
%   run_simple_isac(true,false)       % 仅每个N重采样
%   run_simple_isac(false,true)       % 仅嵌套固定样本
%   run_simple_isac(true,true,true,true) % 双模式+双向干扰均计入

% 设置路径管理
% 获取当前文件所在目录
baseDir = fileparts(mfilename('fullpath'));
% 确保utils目录被添加到路径中，以便调用工具函数
utilsDir = fullfile(baseDir, 'utils');
if ~any(strcmp(path, utilsDir))
    addpath(utilsDir);
end

% 输出目录（按运行时间归档）
outRoot = fullfile(baseDir,'out');
if ~exist(outRoot,'dir'), mkdir(outRoot); end
runStamp = datestr(now,'yyyymmdd_HHMMSS');
outDir = fullfile(outRoot,['run_simple_isac_' runStamp]);
if ~exist(outDir,'dir'), mkdir(outDir); end
fprintf('Results will be saved to: %s\n', outDir);

% 设置随机数种子，保证结果可重复性
rng(1);

% 系统参数配置
M = 6;              % 基站天线数量
K = 4;              % 用户数量
L = 1024;           % 雷达脉冲数
sigmar2 = 10^(-12);  % 雷达接收器噪声功率
sigmak2 = 10^(-12);  % 通信用户噪声功率
sigmat2 = 1;        % 目标反射系数的方差
includeRadarInterferenceInRate = false;  % 设为true可计入雷达对通信的干扰
includeCommInterferenceInRadar = false;  % 设为true可计入通信对雷达的干扰
if numel(varargin) >= 3
    includeRadarInterferenceInRate = logical(varargin{3});
end
if numel(varargin) >= 4
    includeCommInterferenceInRadar = logical(varargin{4});
end
P = 10.^(32/10-3);  % 总发射功率 (32dBm 转换为瓦特)
N_range = 24:8:64;  % RIS单元数量的变化范围
baseline = generate_baseline(M,K);  % 生成基线信道环境
% 新雷达SNR定义包含回程合并增益，按无RIS基线增益缩放门限，保持与旧口径近似可比
gammat_legacy = 10^0.7;
gammat = gammat_legacy*(norm(baseline.hdt)^2);  % 雷达SNR门限值（新口径）

% 双模式开关：
% 1) 每个N重采样（原始行为）
% 2) 嵌套固定样本（先生成max-N，再取前N个单元）
if numel(varargin) >= 1
    enableResampleEachN = logical(varargin{1});
else
    enableResampleEachN = true;
end
if numel(varargin) >= 2
    enableNestedFixedSample = logical(varargin{2});
else
    enableNestedFixedSample = true;
end
Nmax = N_range(end);

if ~(enableResampleEachN || enableNestedFixedSample)
    error('至少需要启用一种模式：enableResampleEachN 或 enableNestedFixedSample。');
end

% 预先生成max-N主样本，用于嵌套固定样本模式
if enableNestedFixedSample
    ris_master = generate_ris_parts(Nmax,baseline);
end

% 初始化性能指标数组（Resample模式）
SR_resample_fixed = zeros(1,length(N_range));
SR_resample_random = zeros(1,length(N_range));
SR_resample_scan_sr = zeros(1,length(N_range));
SR_resample_scan_snr = zeros(1,length(N_range));
SR_resample_noris = zeros(1,length(N_range));

SNR_resample_fixed = zeros(1,length(N_range));
SNR_resample_random = zeros(1,length(N_range));
SNR_resample_scan_sr = zeros(1,length(N_range));
SNR_resample_scan_snr = zeros(1,length(N_range));
SNR_resample_noris = zeros(1,length(N_range));

% 初始化性能指标数组（Nested模式）
SR_nested_fixed = zeros(1,length(N_range));
SR_nested_random = zeros(1,length(N_range));
SR_nested_scan_sr = zeros(1,length(N_range));
SR_nested_scan_snr = zeros(1,length(N_range));
SR_nested_noris = zeros(1,length(N_range));

SNR_nested_fixed = zeros(1,length(N_range));
SNR_nested_random = zeros(1,length(N_range));
SNR_nested_scan_sr = zeros(1,length(N_range));
SNR_nested_scan_snr = zeros(1,length(N_range));
SNR_nested_noris = zeros(1,length(N_range));

% 无RIS基准线与N无关，只需计算一次
Hk_noris = baseline.Hu;
eta_noris = 0;
hrt_noris = 0;
G_noris = zeros(1,size(Hk_noris,2));
phi_noris = 1;
for e = 0:0.05:1
    [Wc_noris_tmp,wr_noris_tmp] = design_w(Hk_noris,baseline.hdt,hrt_noris,G_noris,phi_noris,P,K,e);
    if radar_snr_noris(baseline.hdt,wr_noris_tmp,L,sigmat2,sigmar2,Wc_noris_tmp,includeCommInterferenceInRadar) >= gammat
        eta_noris = e;
        break;
    end
end
[Wc_noris,wr_noris] = design_w(Hk_noris,baseline.hdt,hrt_noris,G_noris,phi_noris,P,K,eta_noris);
SR_noris_val = sum_rate(Hk_noris,Wc_noris,sigmak2,wr_noris,includeRadarInterferenceInRate);
SNR_noris_val = radar_snr_noris(baseline.hdt,wr_noris,L,sigmat2,sigmar2,Wc_noris,includeCommInterferenceInRadar);

% 遍历不同的RIS单元数量，分别评估两种模式
for idx = 1:length(N_range)
    N = N_range(idx);  % 当前RIS单元数量

    if enableResampleEachN
        % 模式A：每个N独立重采样RIS相关信道
        ris = generate_ris_parts(N,baseline);
        Channel.hdt = baseline.hdt;
        Channel.Hu = baseline.Hu;
        Channel.hrt = ris.hrt;
        Channel.G = ris.G;
        Channel.Hru = ris.Hru;

        phi_fixed = compute_phi('fixed_target',Channel);
        phi_random = compute_phi('random',Channel);
        phi_scan_sr = scan_phi('sumrate',Channel,P,K,L,sigmat2,sigmar2,gammat,sigmak2,includeRadarInterferenceInRate,includeCommInterferenceInRadar);
        phi_scan_snr = scan_phi('snr',Channel,P,K,L,sigmat2,sigmar2,gammat,sigmak2,includeRadarInterferenceInRate,includeCommInterferenceInRadar);

        Hk_fixed = Channel.Hu + Channel.Hru*diag(phi_fixed)*Channel.G;
        Hk_random = Channel.Hu + Channel.Hru*diag(phi_random)*Channel.G;
        Hk_scan_sr = Channel.Hu + Channel.Hru*diag(phi_scan_sr)*Channel.G;
        Hk_scan_snr = Channel.Hu + Channel.Hru*diag(phi_scan_snr)*Channel.G;

        eta_fixed = select_eta(gammat,Channel,phi_fixed,P,K,L,sigmat2,sigmar2,includeCommInterferenceInRadar);
        eta_random = select_eta(gammat,Channel,phi_random,P,K,L,sigmat2,sigmar2,includeCommInterferenceInRadar);
        eta_scan_sr = select_eta(gammat,Channel,phi_scan_sr,P,K,L,sigmat2,sigmar2,includeCommInterferenceInRadar);
        eta_scan_snr = select_eta(gammat,Channel,phi_scan_snr,P,K,L,sigmat2,sigmar2,includeCommInterferenceInRadar);

        [Wc_fixed,wr_fixed] = design_w(Hk_fixed,Channel.hdt,Channel.hrt,Channel.G,phi_fixed,P,K,eta_fixed);
        [Wc_random,wr_random] = design_w(Hk_random,Channel.hdt,Channel.hrt,Channel.G,phi_random,P,K,eta_random);
        [Wc_scan_sr,wr_scan_sr] = design_w(Hk_scan_sr,Channel.hdt,Channel.hrt,Channel.G,phi_scan_sr,P,K,eta_scan_sr);
        [Wc_scan_snr,wr_scan_snr] = design_w(Hk_scan_snr,Channel.hdt,Channel.hrt,Channel.G,phi_scan_snr,P,K,eta_scan_snr);

        SR_resample_fixed(idx) = sum_rate(Hk_fixed,Wc_fixed,sigmak2,wr_fixed,includeRadarInterferenceInRate);
        SR_resample_random(idx) = sum_rate(Hk_random,Wc_random,sigmak2,wr_random,includeRadarInterferenceInRate);
        SR_resample_scan_sr(idx) = sum_rate(Hk_scan_sr,Wc_scan_sr,sigmak2,wr_scan_sr,includeRadarInterferenceInRate);
        SR_resample_scan_snr(idx) = sum_rate(Hk_scan_snr,Wc_scan_snr,sigmak2,wr_scan_snr,includeRadarInterferenceInRate);

        SNR_resample_fixed(idx) = radar_snr(Channel.hdt,Channel.hrt,Channel.G,phi_fixed,wr_fixed,L,sigmat2,sigmar2,Wc_fixed,includeCommInterferenceInRadar);
        SNR_resample_random(idx) = radar_snr(Channel.hdt,Channel.hrt,Channel.G,phi_random,wr_random,L,sigmat2,sigmar2,Wc_random,includeCommInterferenceInRadar);
        SNR_resample_scan_sr(idx) = radar_snr(Channel.hdt,Channel.hrt,Channel.G,phi_scan_sr,wr_scan_sr,L,sigmat2,sigmar2,Wc_scan_sr,includeCommInterferenceInRadar);
        SNR_resample_scan_snr(idx) = radar_snr(Channel.hdt,Channel.hrt,Channel.G,phi_scan_snr,wr_scan_snr,L,sigmat2,sigmar2,Wc_scan_snr,includeCommInterferenceInRadar);
    end

    if enableNestedFixedSample
        % 模式B：从max-N主样本切片，形成可比较的嵌套样本
        ris = generate_ris_parts(N,baseline,ris_master);
        Channel.hdt = baseline.hdt;
        Channel.Hu = baseline.Hu;
        Channel.hrt = ris.hrt;
        Channel.G = ris.G;
        Channel.Hru = ris.Hru;

        phi_fixed = compute_phi('fixed_target',Channel);
        phi_random = compute_phi('random',Channel);
        phi_scan_sr = scan_phi('sumrate',Channel,P,K,L,sigmat2,sigmar2,gammat,sigmak2,includeRadarInterferenceInRate,includeCommInterferenceInRadar);
        phi_scan_snr = scan_phi('snr',Channel,P,K,L,sigmat2,sigmar2,gammat,sigmak2,includeRadarInterferenceInRate,includeCommInterferenceInRadar);

        Hk_fixed = Channel.Hu + Channel.Hru*diag(phi_fixed)*Channel.G;
        Hk_random = Channel.Hu + Channel.Hru*diag(phi_random)*Channel.G;
        Hk_scan_sr = Channel.Hu + Channel.Hru*diag(phi_scan_sr)*Channel.G;
        Hk_scan_snr = Channel.Hu + Channel.Hru*diag(phi_scan_snr)*Channel.G;

        eta_fixed = select_eta(gammat,Channel,phi_fixed,P,K,L,sigmat2,sigmar2,includeCommInterferenceInRadar);
        eta_random = select_eta(gammat,Channel,phi_random,P,K,L,sigmat2,sigmar2,includeCommInterferenceInRadar);
        eta_scan_sr = select_eta(gammat,Channel,phi_scan_sr,P,K,L,sigmat2,sigmar2,includeCommInterferenceInRadar);
        eta_scan_snr = select_eta(gammat,Channel,phi_scan_snr,P,K,L,sigmat2,sigmar2,includeCommInterferenceInRadar);

        [Wc_fixed,wr_fixed] = design_w(Hk_fixed,Channel.hdt,Channel.hrt,Channel.G,phi_fixed,P,K,eta_fixed);
        [Wc_random,wr_random] = design_w(Hk_random,Channel.hdt,Channel.hrt,Channel.G,phi_random,P,K,eta_random);
        [Wc_scan_sr,wr_scan_sr] = design_w(Hk_scan_sr,Channel.hdt,Channel.hrt,Channel.G,phi_scan_sr,P,K,eta_scan_sr);
        [Wc_scan_snr,wr_scan_snr] = design_w(Hk_scan_snr,Channel.hdt,Channel.hrt,Channel.G,phi_scan_snr,P,K,eta_scan_snr);

        SR_nested_fixed(idx) = sum_rate(Hk_fixed,Wc_fixed,sigmak2,wr_fixed,includeRadarInterferenceInRate);
        SR_nested_random(idx) = sum_rate(Hk_random,Wc_random,sigmak2,wr_random,includeRadarInterferenceInRate);
        SR_nested_scan_sr(idx) = sum_rate(Hk_scan_sr,Wc_scan_sr,sigmak2,wr_scan_sr,includeRadarInterferenceInRate);
        SR_nested_scan_snr(idx) = sum_rate(Hk_scan_snr,Wc_scan_snr,sigmak2,wr_scan_snr,includeRadarInterferenceInRate);

        SNR_nested_fixed(idx) = radar_snr(Channel.hdt,Channel.hrt,Channel.G,phi_fixed,wr_fixed,L,sigmat2,sigmar2,Wc_fixed,includeCommInterferenceInRadar);
        SNR_nested_random(idx) = radar_snr(Channel.hdt,Channel.hrt,Channel.G,phi_random,wr_random,L,sigmat2,sigmar2,Wc_random,includeCommInterferenceInRadar);
        SNR_nested_scan_sr(idx) = radar_snr(Channel.hdt,Channel.hrt,Channel.G,phi_scan_sr,wr_scan_sr,L,sigmat2,sigmar2,Wc_scan_sr,includeCommInterferenceInRadar);
        SNR_nested_scan_snr(idx) = radar_snr(Channel.hdt,Channel.hrt,Channel.G,phi_scan_snr,wr_scan_snr,L,sigmat2,sigmar2,Wc_scan_snr,includeCommInterferenceInRadar);
    end
end

% 将No-RIS基准线扩展到全N范围（两种模式共用同一条基准线）
SR_resample_noris(:) = SR_noris_val;
SNR_resample_noris(:) = SNR_noris_val;
SR_nested_noris(:) = SR_noris_val;
SNR_nested_noris(:) = SNR_noris_val;

% 绘图：模式A（每个N重采样）
if enableResampleEachN
    fig_sr_resample = figure('Color','w');
    plot(N_range,SR_resample_fixed,'-o','LineWidth',1.5,'Color',[0.8,0,0]); hold on;
    plot(N_range,SR_resample_random,'-s','LineWidth',1.5,'Color',[0,0.5,0]);
    plot(N_range,SR_resample_scan_sr,'-^','LineWidth',1.5,'Color',[0.5,0,0.5]);
    plot(N_range,SR_resample_scan_snr,'-v','LineWidth',1.5,'Color',[0.3,0.3,0.3]);
    plot(N_range,SR_resample_noris,'-d','LineWidth',1.5,'Color',[0,0,0.8]);
    hold off;
    xlabel('N'); ylabel('Sum-rate');
    grid on;
    legend('Fixed','Random','Scan-SR','Scan-SNR','No RIS','Location','best');
    title('Resample each N');

    print(fig_sr_resample, fullfile(outDir,'sumrate_vs_N_resample.png'), '-dpng','-r300');
    print(fig_sr_resample, fullfile(outDir,'sumrate_vs_N_resample.tiff'), '-dtiff','-r300');
    savefig(fig_sr_resample, fullfile(outDir,'sumrate_vs_N_resample.fig'));

    fig_snr_resample = figure('Color','w');
    plot(N_range,10*log10(SNR_resample_fixed),'-o','LineWidth',1.5,'Color',[0.8,0,0]); hold on;
    plot(N_range,10*log10(SNR_resample_random),'-s','LineWidth',1.5,'Color',[0,0.5,0]);
    plot(N_range,10*log10(SNR_resample_scan_sr),'-^','LineWidth',1.5,'Color',[0.5,0,0.5]);
    plot(N_range,10*log10(SNR_resample_scan_snr),'-v','LineWidth',1.5,'Color',[0.3,0.3,0.3]);
    plot(N_range,10*log10(SNR_resample_noris),'-d','LineWidth',1.5,'Color',[0,0,0.8]);
    hold off;
    xlabel('N'); ylabel('Radar SNR (dB)');
    grid on;
    legend('Fixed','Random','Scan-SR','Scan-SNR','No RIS','Location','best');
    title('Resample each N');

    print(fig_snr_resample, fullfile(outDir,'radar_snr_vs_N_resample.png'), '-dpng','-r300');
    print(fig_snr_resample, fullfile(outDir,'radar_snr_vs_N_resample.tiff'), '-dtiff','-r300');
    savefig(fig_snr_resample, fullfile(outDir,'radar_snr_vs_N_resample.fig'));
end

% 绘图：模式B（嵌套固定样本）
if enableNestedFixedSample
    fig_sr_nested = figure('Color','w');
    plot(N_range,SR_nested_fixed,'-o','LineWidth',1.5,'Color',[0.8,0,0]); hold on;
    plot(N_range,SR_nested_random,'-s','LineWidth',1.5,'Color',[0,0.5,0]);
    plot(N_range,SR_nested_scan_sr,'-^','LineWidth',1.5,'Color',[0.5,0,0.5]);
    plot(N_range,SR_nested_scan_snr,'-v','LineWidth',1.5,'Color',[0.3,0.3,0.3]);
    plot(N_range,SR_nested_noris,'-d','LineWidth',1.5,'Color',[0,0,0.8]);
    hold off;
    xlabel('N'); ylabel('Sum-rate');
    grid on;
    legend('Fixed','Random','Scan-SR','Scan-SNR','No RIS','Location','best');
    title('Nested fixed sample');

    print(fig_sr_nested, fullfile(outDir,'sumrate_vs_N_nested.png'), '-dpng','-r300');
    print(fig_sr_nested, fullfile(outDir,'sumrate_vs_N_nested.tiff'), '-dtiff','-r300');
    savefig(fig_sr_nested, fullfile(outDir,'sumrate_vs_N_nested.fig'));

    fig_snr_nested = figure('Color','w');
    plot(N_range,10*log10(SNR_nested_fixed),'-o','LineWidth',1.5,'Color',[0.8,0,0]); hold on;
    plot(N_range,10*log10(SNR_nested_random),'-s','LineWidth',1.5,'Color',[0,0.5,0]);
    plot(N_range,10*log10(SNR_nested_scan_sr),'-^','LineWidth',1.5,'Color',[0.5,0,0.5]);
    plot(N_range,10*log10(SNR_nested_scan_snr),'-v','LineWidth',1.5,'Color',[0.3,0.3,0.3]);
    plot(N_range,10*log10(SNR_nested_noris),'-d','LineWidth',1.5,'Color',[0,0,0.8]);
    hold off;
    xlabel('N'); ylabel('Radar SNR (dB)');
    grid on;
    legend('Fixed','Random','Scan-SR','Scan-SNR','No RIS','Location','best');
    title('Nested fixed sample');

    print(fig_snr_nested, fullfile(outDir,'radar_snr_vs_N_nested.png'), '-dpng','-r300');
    print(fig_snr_nested, fullfile(outDir,'radar_snr_vs_N_nested.tiff'), '-dtiff','-r300');
    savefig(fig_snr_nested, fullfile(outDir,'radar_snr_vs_N_nested.fig'));
end

% 额外对比图：仅比较Fixed策略在两种模式下的差异
if enableResampleEachN && enableNestedFixedSample
    fig_fixed_sr_compare = figure('Color','w');
    plot(N_range,SR_resample_fixed,'-o','LineWidth',1.5,'Color',[0.8,0,0]); hold on;
    plot(N_range,SR_nested_fixed,'-s','LineWidth',1.5,'Color',[0,0.5,0]);
    plot(N_range,SR_resample_noris,'-d','LineWidth',1.5,'Color',[0,0,0.8]);
    hold off;
    xlabel('N'); ylabel('Sum-rate');
    grid on;
    legend('Fixed (Resample each N)','Fixed (Nested sample)','No RIS','Location','best');
    title('Fixed phase: channel sampling mode comparison');

    print(fig_fixed_sr_compare, fullfile(outDir,'fixed_sumrate_mode_compare.png'), '-dpng','-r300');
    print(fig_fixed_sr_compare, fullfile(outDir,'fixed_sumrate_mode_compare.tiff'), '-dtiff','-r300');
    savefig(fig_fixed_sr_compare, fullfile(outDir,'fixed_sumrate_mode_compare.fig'));

    fig_fixed_snr_compare = figure('Color','w');
    plot(N_range,10*log10(SNR_resample_fixed),'-o','LineWidth',1.5,'Color',[0.8,0,0]); hold on;
    plot(N_range,10*log10(SNR_nested_fixed),'-s','LineWidth',1.5,'Color',[0,0.5,0]);
    plot(N_range,10*log10(SNR_resample_noris),'-d','LineWidth',1.5,'Color',[0,0,0.8]);
    hold off;
    xlabel('N'); ylabel('Radar SNR (dB)');
    grid on;
    legend('Fixed (Resample each N)','Fixed (Nested sample)','No RIS','Location','best');
    title('Fixed phase: channel sampling mode comparison');

    print(fig_fixed_snr_compare, fullfile(outDir,'fixed_radarsnr_mode_compare.png'), '-dpng','-r300');
    print(fig_fixed_snr_compare, fullfile(outDir,'fixed_radarsnr_mode_compare.tiff'), '-dtiff','-r300');
    savefig(fig_fixed_snr_compare, fullfile(outDir,'fixed_radarsnr_mode_compare.fig'));
end

% 保存关键数组，便于后处理和论文作图
save(fullfile(outDir,'mode_compare_results.mat'), ...
    'N_range', ...
    'SR_resample_fixed','SR_resample_random','SR_resample_scan_sr','SR_resample_scan_snr','SR_resample_noris', ...
    'SNR_resample_fixed','SNR_resample_random','SNR_resample_scan_sr','SNR_resample_scan_snr','SNR_resample_noris', ...
    'SR_nested_fixed','SR_nested_random','SR_nested_scan_sr','SR_nested_scan_snr','SR_nested_noris', ...
    'SNR_nested_fixed','SNR_nested_random','SNR_nested_scan_sr','SNR_nested_scan_snr','SNR_nested_noris', ...
    'enableResampleEachN','enableNestedFixedSample');

% 全局相位变化对性能的影响分析部分
N0 = N_range(end);  % 使用最大的RIS单元数量
if enableNestedFixedSample
    ris0 = generate_ris_parts(N0,baseline,ris_master);
else
    ris0 = generate_ris_parts(N0,baseline);  % 基于基线生成RIS部分
end
Channel0.hdt = baseline.hdt; Channel0.Hu = baseline.Hu;
Channel0.hrt = ris0.hrt; Channel0.G = ris0.G; Channel0.Hru = ris0.Hru;
phi_base = compute_phi('fixed_target',Channel0);  % 基于目标的基础相位

% 生成多个全局相位偏移点
rhdeltas = linspace(0,2*pi,30);  % 0到2π的30个均匀分布点
SR_step = zeros(1,length(rhdeltas));  % 存储不同相位下的速率
SNR_step = zeros(1,length(rhdeltas)); % 存储不同相位下的SNR

% 计算每个全局相位偏移下的系统性能
for t = 1:length(rhdeltas)
    % 应用全局相位偏移
    phi_step = phi_base.*exp(1i*rhdeltas(t));
    % 计算等效信道
    Hk_step = Channel0.Hu + Channel0.Hru*diag(phi_step)*Channel0.G;
    % 选择功率分配系数
    eta_step = select_eta(gammat,Channel0,phi_step,P,K,L,sigmat2,sigmar2,includeCommInterferenceInRadar);
    % 设计波束赋形向量
    [Wc_step,wr_step] = design_w(Hk_step,Channel0.hdt,Channel0.hrt,Channel0.G,phi_step,P,K,eta_step);
    % 计算性能指标
    SR_step(t) = sum_rate(Hk_step,Wc_step,sigmak2,wr_step,includeRadarInterferenceInRate);
    SNR_step(t) = radar_snr(Channel0.hdt,Channel0.hrt,Channel0.G,phi_step,wr_step,L,sigmat2,sigmar2,Wc_step,includeCommInterferenceInRadar);
end

% 绘制全局相位对通信速率的影响并保存
fig3 = figure('Color','w'); 
plot(rhdeltas,SR_step,'-o','LineWidth',1.5,'Color',[0.8,0,0]); 
grid on;
xlabel('Global phase (rad)'); ylabel('Sum-rate');
% 设置横坐标刻度为0-2π，间隔π/4，并使用数学符号标签
set(gca,'XTick', 0:pi/4:2*pi);
set(gca,'XTickLabel', {'0', 'π/4', 'π/2', '3π/4', 'π', '5π/4', '3π/2', '7π/4', '2π'});

print(fig3, fullfile(outDir,'global_phase_sumrate.png'), '-dpng','-r300');
print(fig3, fullfile(outDir,'global_phase_sumrate.tiff'), '-dtiff','-r300');
savefig(fig3, fullfile(outDir,'global_phase_sumrate.fig'));

% 绘制全局相位对雷达SNR的影响并保存
fig4 = figure('Color','w'); 
plot(rhdeltas,10*log10(SNR_step),'-o','LineWidth',1.5,'Color',[0,0.5,0]); 
grid on;
xlabel('Global phase (rad)'); ylabel('Radar SNR (dB)');
% 设置横坐标刻度
set(gca,'XTick', 0:pi/4:2*pi);
set(gca,'XTickLabel', {'0', 'π/4', 'π/2', '3π/4', 'π', '5π/4', '3π/2', '7π/4', '2π'});

print(fig4, fullfile(outDir,'global_phase_radarsnr.png'), '-dpng','-r300');
print(fig4, fullfile(outDir,'global_phase_radarsnr.tiff'), '-dtiff','-r300');
savefig(fig4, fullfile(outDir,'global_phase_radarsnr.fig'));

% η权衡图（上下堆放）- 分析功率分配系数对通信和雷达性能的影响
N0 = N_range(end);
if enableNestedFixedSample
    ris0 = generate_ris_parts(N0,baseline,ris_master);
else
    ris0 = generate_ris_parts(N0,baseline);
end
Channel0.hdt = baseline.hdt; Channel0.Hu = baseline.Hu;
Channel0.hrt = ris0.hrt; Channel0.G = ris0.G; Channel0.Hru = ris0.Hru;
phi = compute_phi('fixed_target',Channel0);
Hk = Channel0.Hu + Channel0.Hru*diag(phi)*Channel0.G;
etas = 0:0.05:1;  % 功率分配系数从0到1，步长0.05
SR_eta = zeros(size(etas)); SNR_eta = zeros(size(etas));

% 计算不同η值下的系统性能
for i = 1:numel(etas)
    e = etas(i);
    [Wc_eta,wr_eta] = design_w(Hk,Channel0.hdt,Channel0.hrt,Channel0.G,phi,P,K,e);
    SR_eta(i) = sum_rate(Hk,Wc_eta,sigmak2,wr_eta,includeRadarInterferenceInRate);
    SNR_eta(i) = radar_snr(Channel0.hdt,Channel0.hrt,Channel0.G,phi,wr_eta,L,sigmat2,sigmar2,Wc_eta,includeCommInterferenceInRadar);
end

% 创建η权衡图，设置图形属性和位置
fig_eta = figure('Color','w');
set(fig_eta,'Position',[120 120 640 800]);  % 设置图形窗口位置和大小

% 上子图：雷达SNR与η的关系
subplot(2,1,1);
plot(etas,10*log10(SNR_eta),'-','LineWidth',1.5,'Color',[0.8,0,0]); hold on; grid on;
% 标记关键η取值点（0.2, 0.5, 0.8）
[~,i1] = min(abs(etas-0.2)); [~,i2] = min(abs(etas-0.5)); [~,i3] = min(abs(etas-0.8));
plot(etas(i1),10*log10(SNR_eta(i1)),'o','MarkerSize',7,'MarkerFaceColor',[0.8,0,0],'MarkerEdgeColor',[0.8,0,0]);
plot(etas(i2),10*log10(SNR_eta(i2)),'s','MarkerSize',7,'MarkerFaceColor',[0,0.5,0],'MarkerEdgeColor',[0,0,0.5]);
plot(etas(i3),10*log10(SNR_eta(i3)),'^','MarkerSize',7,'MarkerFaceColor',[0,0,0.8],'MarkerEdgeColor',[0,0,0.8]);
xlabel('\eta'); ylabel('Radar SNR (dB)');
legend('\eta=0.2','\eta=0.5','\eta=0.8','Location','best');
xlim([0 1]);

% 标注通信优先和感知优先区域
yl = ylim; patch([0 0.4 0.4 0],[yl(1) yl(1) yl(2) yl(2)], [0.9 0.9 0.9],'FaceAlpha',0.2,'EdgeColor','none');
text(0.2,yl(1)+(yl(2)-yl(1))*0.1,'通信优先','HorizontalAlignment','center');
patch([0.6 1 1 0.6],[yl(1) yl(1) yl(2) yl(2)], [0.9 0.9 0.9],'FaceAlpha',0.2,'EdgeColor','none');
text(0.8,yl(1)+(yl(2)-yl(1))*0.1,'感知优先','HorizontalAlignment','center');

% 下子图：通信速率与η的关系
subplot(2,1,2);
plot(etas,SR_eta,'-','LineWidth',1.5,'Color',[0,0.5,0]); hold on; grid on;
% 标记相同的关键η取值点
plot(etas(i1),SR_eta(i1),'o','MarkerSize',7,'MarkerFaceColor',[0.8,0,0],'MarkerEdgeColor',[0.8,0,0]);
plot(etas(i2),SR_eta(i2),'s','MarkerSize',7,'MarkerFaceColor',[0,0,0.5],'MarkerEdgeColor',[0,0,0.5]);
plot(etas(i3),SR_eta(i3),'^','MarkerSize',7,'MarkerFaceColor',[0,0,0.8],'MarkerEdgeColor',[0,0,0.8]);
xlabel('\eta'); ylabel('Sum-rate (bps/Hz)');
legend('\eta=0.2','\eta=0.5','\eta=0.8','Location','best');
xlim([0 1]);

% 同样标注通信优先和感知优先区域
yl2 = ylim; patch([0 0.4 0.4 0],[yl2(1) yl2(1) yl2(2) yl2(2)], [0.9 0.9 0.9],'FaceAlpha',0.2,'EdgeColor','none');
text(0.2,yl2(1)+(yl2(2)-yl2(1))*0.1,'通信优先','HorizontalAlignment','center');
patch([0.6 1 1 0.6],[yl2(1) yl2(1) yl2(2) yl2(2)], [0.9 0.9 0.9],'FaceAlpha',0.2,'EdgeColor','none');
text(0.8,yl2(1)+(yl2(2)-yl2(1))*0.1,'感知优先','HorizontalAlignment','center');

% 保存结果图片
print(fig_eta, fullfile(outDir,'eta_tradeoff.png'), '-dpng','-r300');
print(fig_eta, fullfile(outDir,'eta_tradeoff.tiff'), '-dtiff','-r300');
savefig(fig_eta, fullfile(outDir,'eta_tradeoff.fig'));
end