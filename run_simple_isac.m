function run_simple_isac()
% 主函数：运行RIS辅助的ISAC(集成感知与通信)系统仿真
% 功能：对比不同RIS相位配置策略下的通信速率和雷达性能

% 设置路径管理
% 获取当前文件所在目录
baseDir = fileparts(mfilename('fullpath'));
% 确保utils目录被添加到路径中，以便调用工具函数
utilsDir = fullfile(baseDir, 'utils');
if ~any(strcmp(path, utilsDir))
    addpath(utilsDir);
end

% 输出目录（保存所有图像）
outDir = fullfile(baseDir,'out');
if ~exist(outDir,'dir'), mkdir(outDir); end

% 设置随机数种子，保证结果可重复性
rng(1);

% 系统参数配置
M = 6;              % 基站天线数量
K = 4;              % 用户数量
L = 1024;           % 雷达脉冲数
sigmar2 = 10^(-12);  % 雷达接收器噪声功率
sigmak2 = 10^(-12);  % 通信用户噪声功率
sigmat2 = 1;        % 目标反射系数的方差
P = 10.^(32/10-3);  % 总发射功率 (32dBm 转换为瓦特)
N_range = 24:8:64;  % RIS单元数量的变化范围
gammat = 10^0.7;  % 雷达SNR门限值
baseline = generate_baseline(M,K);  % 生成基线信道环境

% 初始化性能指标数组 - 用于存储不同RIS相位策略下的性能
SR_fixed = zeros(1,length(N_range));    % 固定相位策略的总和速率
SR_random = zeros(1,length(N_range));   % 随机相位策略的总和速率
SR_scan_sr = zeros(1,length(N_range));  % 扫描优化速率策略的总和速率
SR_scan_snr = zeros(1,length(N_range)); % 扫描优化SNR策略的总和速率
SR_noris = zeros(1,length(N_range));    % 无RIS场景的总和速率

SNR_fixed = zeros(1,length(N_range));   % 固定相位策略的雷达SNR
SNR_random = zeros(1,length(N_range));  % 随机相位策略的雷达SNR
SNR_scan_sr = zeros(1,length(N_range)); % 扫描优化速率策略的雷达SNR
SNR_scan_snr = zeros(1,length(N_range));% 扫描优化SNR策略的雷达SNR
SNR_noris = zeros(1,length(N_range));   % 无RIS场景的雷达SNR

% 遍历不同的RIS单元数量，评估系统性能随RIS数量的变化
for idx = 1:length(N_range)
    N = N_range(idx);  % 当前RIS单元数量
    
    % 生成RIS相关部分，并与固定基线组合
    ris = generate_ris_parts(N,baseline);
    Channel.hdt = baseline.hdt;  % 基站-目标直接信道
    Channel.Hu = baseline.Hu;    % 基站-用户直接信道
    Channel.hrt = ris.hrt;       % RIS-目标信道
    Channel.G = ris.G;           % 基站-RIS信道
    Channel.Hru = ris.Hru;       % RIS-用户信道
    
    % 计算不同相位策略下的RIS相位矩阵
    phi_fixed = compute_phi('fixed_target',Channel);    % 固定目标相位
    phi_random = compute_phi('random',Channel);         % 随机相位
    phi_scan_sr = scan_phi('sumrate',Channel,P,K,L,sigmat2,sigmar2,gammat,sigmak2); % 优化速率的扫描相位
    phi_scan_snr = scan_phi('snr',Channel,P,K,L,sigmat2,sigmar2,gammat,sigmak2);    % 优化SNR的扫描相位
    
    % 计算等效信道矩阵（考虑RIS的影响）
    Hk_fixed = Channel.Hu + Channel.Hru*diag(phi_fixed)*Channel.G;
    Hk_random = Channel.Hu + Channel.Hru*diag(phi_random)*Channel.G;
    Hk_scan_sr = Channel.Hu + Channel.Hru*diag(phi_scan_sr)*Channel.G;
    Hk_scan_snr = Channel.Hu + Channel.Hru*diag(phi_scan_snr)*Channel.G;
    
    % 选择功率分配系数eta（雷达与通信之间的功率分配）
    eta_fixed = select_eta(gammat,Channel,phi_fixed,P,K,L,sigmat2,sigmar2);
    eta_random = select_eta(gammat,Channel,phi_random,P,K,L,sigmat2,sigmar2);
    eta_scan_sr = select_eta(gammat,Channel,phi_scan_sr,P,K,L,sigmat2,sigmar2);
    eta_scan_snr = select_eta(gammat,Channel,phi_scan_snr,P,K,L,sigmat2,sigmar2);
    
    % 设计波束赋形向量（通信和雷达）
    [Wc_fixed,wr_fixed] = design_w(Hk_fixed,Channel.hdt,P,K,eta_fixed);
    [Wc_random,wr_random] = design_w(Hk_random,Channel.hdt,P,K,eta_random);
    [Wc_scan_sr,wr_scan_sr] = design_w(Hk_scan_sr,Channel.hdt,P,K,eta_scan_sr);
    [Wc_scan_snr,wr_scan_snr] = design_w(Hk_scan_snr,Channel.hdt,P,K,eta_scan_snr);
    
    % 计算通信性能指标（总和速率）
    SR_fixed(idx) = sum_rate(Hk_fixed,Wc_fixed,sigmak2);
    SR_random(idx) = sum_rate(Hk_random,Wc_random,sigmak2);
    SR_scan_sr(idx) = sum_rate(Hk_scan_sr,Wc_scan_sr,sigmak2);
    SR_scan_snr(idx) = sum_rate(Hk_scan_snr,Wc_scan_snr,sigmak2);
    
    % 计算雷达性能指标（SNR）
    SNR_fixed(idx) = radar_snr(Channel.hdt,Channel.hrt,Channel.G,phi_fixed,wr_fixed,L,sigmat2,sigmar2);
    SNR_random(idx) = radar_snr(Channel.hdt,Channel.hrt,Channel.G,phi_random,wr_random,L,sigmat2,sigmar2);
    SNR_scan_sr(idx) = radar_snr(Channel.hdt,Channel.hrt,Channel.G,phi_scan_sr,wr_scan_sr,L,sigmat2,sigmar2);
    SNR_scan_snr(idx) = radar_snr(Channel.hdt,Channel.hrt,Channel.G,phi_scan_snr,wr_scan_snr,L,sigmat2,sigmar2);
    
    % 无RIS场景的性能计算（使用固定基线环境，独立于N）
    Hk_noris = baseline.Hu;  % 无RIS时只有直接信道
    eta_noris = 0;           % 初始化功率分配系数
    % 搜索满足雷达SNR门限的最小通信功率分配比例
    for e = 0:0.05:0.5
        [~,wr_noris_tmp] = design_w(Hk_noris,baseline.hdt,P,K,e);
        if radar_snr_noris(baseline.hdt,wr_noris_tmp,L,sigmat2,sigmar2) >= gammat
            eta_noris = e; break;
        end
    end
    [Wc_noris,wr_noris] = design_w(Hk_noris,baseline.hdt,P,K,eta_noris);
    SR_noris(idx) = sum_rate(Hk_noris,Wc_noris,sigmak2);
    SNR_noris(idx) = radar_snr_noris(baseline.hdt,wr_noris,L,sigmat2,sigmar2);
end

% 绘制通信速率随RIS单元数量变化的曲线并保存
fig1 = figure('Color','w'); 
plot(N_range,SR_fixed,'-o','LineWidth',1.5,'Color',[0.8,0,0]); hold on;
plot(N_range,SR_random,'-s','LineWidth',1.5,'Color',[0,0.5,0]);
plot(N_range,SR_scan_sr,'-^','LineWidth',1.5,'Color',[0.5,0,0.5]);
plot(N_range,SR_scan_snr,'-v','LineWidth',1.5,'Color',[0.3,0.3,0.3]);
plot(N_range,SR_noris,'-d','LineWidth',1.5,'Color',[0,0,0.8]); 
hold off;
xlabel('N'); ylabel('Sum-rate'); 
grid on; 
legend('Fixed','Random','Scan-SR','Scan-SNR','No RIS');

print(fig1, fullfile(outDir,'sumrate_vs_N.png'), '-dpng','-r300');
print(fig1, fullfile(outDir,'sumrate_vs_N.tiff'), '-dtiff','-r300');
savefig(fig1, fullfile(outDir,'sumrate_vs_N.fig'));

% 绘制雷达SNR随RIS单元数量变化的曲线（dB表示）并保存
fig2 = figure('Color','w'); 
plot(N_range,10*log10(SNR_fixed),'-o','LineWidth',1.5,'Color',[0.8,0,0]); hold on;
plot(N_range,10*log10(SNR_random),'-s','LineWidth',1.5,'Color',[0,0.5,0]);
plot(N_range,10*log10(SNR_scan_sr),'-^','LineWidth',1.5,'Color',[0.5,0,0.5]);
plot(N_range,10*log10(SNR_scan_snr),'-v','LineWidth',1.5,'Color',[0.3,0.3,0.3]);
plot(N_range,10*log10(SNR_noris),'-d','LineWidth',1.5,'Color',[0,0,0.8]); 
hold off;
xlabel('N'); ylabel('Radar SNR (dB)'); 
grid on; 
legend('Fixed','Random','Scan-SR','Scan-SNR','No RIS');

print(fig2, fullfile(outDir,'radar_snr_vs_N.png'), '-dpng','-r300');
print(fig2, fullfile(outDir,'radar_snr_vs_N.tiff'), '-dtiff','-r300');
savefig(fig2, fullfile(outDir,'radar_snr_vs_N.fig'));

% 全局相位变化对性能的影响分析部分
N0 = N_range(end);  % 使用最大的RIS单元数量
ris0 = generate_ris_parts(N0,baseline);  % 基于基线生成RIS部分
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
    eta_step = select_eta(gammat,Channel0,phi_step,P,K,L,sigmat2,sigmar2);
    % 设计波束赋形向量
    [Wc_step,wr_step] = design_w(Hk_step,Channel0.hdt,P,K,eta_step);
    % 计算性能指标
    SR_step(t) = sum_rate(Hk_step,Wc_step,sigmak2);
    SNR_step(t) = radar_snr(Channel0.hdt,Channel0.hrt,Channel0.G,phi_step,wr_step,L,sigmat2,sigmar2);
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
ris0 = generate_ris_parts(N0,baseline);
Channel0.hdt = baseline.hdt; Channel0.Hu = baseline.Hu;
Channel0.hrt = ris0.hrt; Channel0.G = ris0.G; Channel0.Hru = ris0.Hru;
phi = compute_phi('fixed_target',Channel0);
Hk = Channel0.Hu + Channel0.Hru*diag(phi)*Channel0.G;
etas = 0:0.05:1;  % 功率分配系数从0到1，步长0.05
SR_eta = zeros(size(etas)); SNR_eta = zeros(size(etas));

% 计算不同η值下的系统性能
for i = 1:numel(etas)
    e = etas(i);
    [Wc_eta,wr_eta] = design_w(Hk,Channel0.hdt,P,K,e);
    SR_eta(i) = sum_rate(Hk,Wc_eta,sigmak2);
    SNR_eta(i) = radar_snr(Channel0.hdt,Channel0.hrt,Channel0.G,phi,wr_eta,L,sigmat2,sigmar2);
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