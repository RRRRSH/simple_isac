function plot_eta_tradeoff()
% plot_eta_tradeoff - 绘制功率分配参数η对ISAC系统性能的影响
% 
% 该函数生成通信-雷达功率分配参数η的权衡关系图，展示不同η值下通信速率和雷达SNR的变化
% 固定RIS单元数N=64，采用fixed-target相位策略
% 结果将保存为PNG和TIFF格式图像

% 设置工作环境
baseDir = fileparts(mfilename('fullpath'));  % 获取当前脚本所在目录
addpath(fullfile(baseDir,'utils'));  % 添加工具函数目录到路径
rng(1);  % 设置随机数种子，确保结果可重复

% 输出目录（按运行时间归档）
outRoot = fullfile(baseDir,'out');
if ~exist(outRoot,'dir'), mkdir(outRoot); end
runStamp = datestr(now,'yyyymmdd_HHMMSS');
outdir = fullfile(outRoot,['plot_eta_tradeoff_' runStamp]);
if ~exist(outdir,'dir'), mkdir(outdir); end
fprintf('Results will be saved to: %s\n', outdir);

% 系统参数配置
M = 6;  % 基站天线数
K = 4;  % 用户数量
L = 1024;  % 雷达脉冲数
sigmar2 = 10^(-12);  % 雷达接收噪声功率
sigmak2 = 10^(-12);  % 通信用户噪声功率
sigmat2 = 1;  % 目标反射系数方差
includeRadarInterferenceInRate = false;  % 设为true可计入雷达对通信的干扰
includeCommInterferenceInRadar = false;  % 设为true可计入通信对雷达的干扰
P = 10.^(32/10-3);  % 总发射功率 (32dBm转换为瓦特)
N0 = 64;  % RIS单元数

% 生成系统模型参数
baseline = generate_baseline(M,K);  % 生成无RIS情况下的基线信道参数
ris0 = generate_ris_parts(N0,baseline);  % 生成RIS相关的信道参数

% 组装完整的信道参数结构
Channel0.hdt = baseline.hdt;  % 基站-目标直接信道 (M×1)
Channel0.Hu = baseline.Hu;     % 基站-用户直接信道 (K×M)
Channel0.hrt = ris0.hrt;       % RIS-目标信道 (N×1)
Channel0.G = ris0.G;           % 用户-RIS信道 (N×K)
Channel0.Hru = ris0.Hru;       % 基站-RIS信道 (M×N)

% 计算RIS相位矩阵
phi = compute_phi('fixed_target',Channel0);  % 使用基于目标的固定相位策略

% 计算等效信道矩阵 (考虑RIS反射)
Hk = Channel0.Hu + Channel0.Hru*diag(phi)*Channel0.G;

% 生成不同的功率分配参数η值 (0到1，步长0.05)
etas = 0:0.05:1;  % η表示分配给雷达的功率比例

% 初始化存储性能指标的数组
SR = zeros(size(etas));  % 存储各η值对应的通信总速率
SNR = zeros(size(etas));  % 存储各η值对应的雷达SNR

% 遍历所有η值，计算对应性能
for i = 1:numel(etas)
    e = etas(i);  % 当前功率分配参数
    % 设计波束赋形向量
    [Wc,wr] = design_w(Hk,Channel0.hdt,Channel0.hrt,Channel0.G,phi,P,K,e);  % Wc:通信波束, wr:雷达波束
    % 计算通信性能指标
    SR(i) = sum_rate(Hk,Wc,sigmak2,wr,includeRadarInterferenceInRate);  % 计算通信总速率
    % 计算雷达性能指标
    SNR(i) = radar_snr(Channel0.hdt,Channel0.hrt,Channel0.G,phi,wr,L,sigmat2,sigmar2,Wc,includeCommInterferenceInRadar);
end

% 创建图形并设置属性
fig = figure('Color','w');  % 创建白色背景的图形窗口
set(fig,'Position',[100 100 640 800]);  % 设置图形窗口位置和大小

% 绘制雷达SNR随η变化的曲线（上图）
subplot(2,1,1);  % 创建2行1列的子图，当前为第1个子图
plot(etas,10*log10(SNR),'-','LineWidth',1.5,'Color',[0.8,0,0]);  % 绘制雷达SNR曲线(dB)
hold on; grid on;  % 保持当前图形并添加网格

% 标记特定η值对应的点(0.2, 0.5, 0.8)
[~,i1] = min(abs(etas-0.2));  % 找到最接近η=0.2的索引
[~,i2] = min(abs(etas-0.5));  % 找到最接近η=0.5的索引
[~,i3] = min(abs(etas-0.8));  % 找到最接近η=0.8的索引

% 绘制特定点并设置不同标记样式
plot(etas(i1),10*log10(SNR(i1)),'o','MarkerSize',7,'MarkerFaceColor',[0.8,0,0],'MarkerEdgeColor',[0.8,0,0]);
plot(etas(i2),10*log10(SNR(i2)),'s','MarkerSize',7,'MarkerFaceColor',[0,0.5,0],'MarkerEdgeColor',[0,0.5,0]);
plot(etas(i3),10*log10(SNR(i3)),'^','MarkerSize',7,'MarkerFaceColor',[0,0,0.8],'MarkerEdgeColor',[0,0,0.8]);

% 设置坐标轴标签和图例
xlabel('\eta');  % x轴：功率分配参数η
ylabel('Radar SNR (dB)');  % y轴：雷达信噪比(dB)
legend('\eta=0.2','\eta=0.5','\eta=0.8','Location','best');  % 添加图例
xlim([0 1]);  % 设置x轴范围

% 添加功能区标注（通信优先和感知优先区域）
yl = ylim;  % 获取当前y轴范围
% 通信优先区域阴影（η较小）
patch([0 0.4 0.4 0],[yl(1) yl(1) yl(2) yl(2)], [0.9 0.9 0.9],'FaceAlpha',0.2,'EdgeColor','none');
text(0.2,yl(1)+(yl(2)-yl(1))*0.1,'通信优先','HorizontalAlignment','center');  % 添加文字标签
% 感知优先区域阴影（η较大）
patch([0.6 1 1 0.6],[yl(1) yl(1) yl(2) yl(2)], [0.9 0.9 0.9],'FaceAlpha',0.2,'EdgeColor','none');
text(0.8,yl(1)+(yl(2)-yl(1))*0.1,'感知优先','HorizontalAlignment','center');  % 添加文字标签

% 绘制通信速率随η变化的曲线（下图）
subplot(2,1,2);  % 创建2行1列的子图，当前为第2个子图
plot(etas,SR,'-','LineWidth',1.5,'Color',[0,0.5,0]);  % 绘制通信总速率曲线
hold on; grid on;  % 保持当前图形并添加网格

% 在通信速率图上标记相同的特定η值点
plot(etas(i1),SR(i1),'o','MarkerSize',7,'MarkerFaceColor',[0.8,0,0],'MarkerEdgeColor',[0.8,0,0]);
plot(etas(i2),SR(i2),'s','MarkerSize',7,'MarkerFaceColor',[0,0.5,0],'MarkerEdgeColor',[0,0.5,0]);
plot(etas(i3),SR(i3),'^','MarkerSize',7,'MarkerFaceColor',[0,0,0.8],'MarkerEdgeColor',[0,0,0.8]);

% 设置坐标轴标签和图例
xlabel('\eta');  % x轴：功率分配参数η
ylabel('Sum-rate (bps/Hz)');  % y轴：通信总速率(bps/Hz)
legend('\eta=0.2','\eta=0.5','\eta=0.8','Location','best');  % 添加图例
xlim([0 1]);  % 设置x轴范围

% 添加功能区标注
yl2 = ylim;  % 获取当前y轴范围
% 通信优先区域阴影
patch([0 0.4 0.4 0],[yl2(1) yl2(1) yl2(2) yl2(2)], [0.9 0.9 0.9],'FaceAlpha',0.2,'EdgeColor','none');
text(0.2,yl2(1)+(yl2(2)-yl2(1))*0.1,'通信优先','HorizontalAlignment','center');
% 感知优先区域阴影
patch([0.6 1 1 0.6],[yl2(1) yl2(1) yl2(2) yl2(2)], [0.9 0.9 0.9],'FaceAlpha',0.2,'EdgeColor','none');
text(0.8,yl2(1)+(yl2(2)-yl2(1))*0.1,'感知优先','HorizontalAlignment','center');

% 添加图形说明文本框
annotation('textbox',[0.12 0.01 0.76 0.08],'String', ...
    '固定N=64，采用fixed-target相位；η为雷达功率比例。低η通信优先（左侧阴影），高η感知优先（右侧阴影）。', ...
    'HorizontalAlignment','center','EdgeColor','none');  % 居中显示，无边框

% 保存图形到输出目录
% 保存为PNG格式，分辨率300dpi
print(fig, fullfile(outdir,'eta_tradeoff.png'), '-dpng','-r300');
% 保存为TIFF格式，分辨率300dpi
print(fig, fullfile(outdir,'eta_tradeoff.tiff'), '-dtiff','-r300');

end