function run_noris_full()
% run_noris_full
% 仅运行 No-RIS 的完整流程：
% 1) 固定发射功率 P=1W
% 2) 在雷达SNR门限约束下优化 No-RIS 波束
% 3) 生成并保存 No-RIS 对应图像与数据

baseDir = fileparts(mfilename('fullpath'));
if isempty(baseDir), baseDir = pwd; end
addpath(fullfile(baseDir,'utils'));

fprintf('\n========== No-RIS Full Run ==========%s\n', newline);

% 输出目录（按运行时间归档）
outRoot = fullfile(baseDir,'out');
if ~exist(outRoot,'dir'), mkdir(outRoot); end
runStamp = datestr(now,'yyyymmdd_HHMMSS');
outDir = fullfile(outRoot,['run_noris_full_' runStamp]);
if ~exist(outDir,'dir'), mkdir(outDir); end
fprintf('Results will be saved to: %s\n', outDir);

fprintf('[1/8] 初始化参数...\n');
rng(1);
M = 6;
K = 4;
L = 1024;
sigmar2 = 1e-12;
sigmak2 = 1e-12;
sigmat2 = 1;
P = 1.0;  % 固定1W

% 统一雷达SNR门限：固定为7 dB（线性值）
gammat = 10^(7/10);
fprintf('      M=%d, K=%d, L=%d, P=%.3f W\n', M, K, L, P);
fprintf('      gammat(linear)=%.3e, gammat(dB)=%.3f dB\n', gammat, 10*log10(gammat));

fprintf('[2/8] 生成基线信道...\n');
baseline = generate_baseline(M,K);

fprintf('[3/8] 构造No-RIS信道并优化波束（输出迭代进度）...\n');
Channel.hdt = baseline.hdt;
Channel.Hu = baseline.Hu;
Channel.hrt = 0;
Channel.G = zeros(1,M);
Channel.Hru = zeros(K,1);
phi_noris = 1;

[Wc_opt,wr_opt,~,Vsr_opt,gammat_opt] = optimize_w_for_fixed_phi( ...
    Channel,phi_noris,P,K,L,sigmat2,sigmar2,sigmak2,gammat,true);

fprintf('[4/8] 计算核心性能指标...\n');
Hk_noris = Channel.Hu;
SR = sum_rate(Hk_noris,Wc_opt,sigmak2,wr_opt,false);
SNR_legacy = radar_snr_noris(Channel.hdt,wr_opt,L,sigmat2,sigmar2,Wc_opt,false);
ptotal = norm(Wc_opt,'fro')^2 + norm(wr_opt,2)^2;

fprintf('[5/8] 生成No-RIS空间功率热力图...\n');
scene = build_paper_scene(K);
allPts = [scene.bs; scene.ris; scene.target; scene.users];
margin = 10;
x_range = floor(min(allPts(:,1))-margin):0.5:ceil(max(allPts(:,1))+margin);
y_range = floor(min(allPts(:,2))-margin):0.5:ceil(max(allPts(:,2))+margin);
[X,Y] = meshgrid(x_range,y_range);

alpha_bt = 2.2;
dmin_bs = 1e-6;

P_rad_NoRIS = zeros(size(X));
wr_radar_sum = sum(wr_opt,2);
w_total_noris = wr_radar_sum + sum(Wc_opt,2);
P_total_NoRIS = zeros(size(X));

for i = 1:size(X,1)
    for j = 1:size(X,2)
        pos = [X(i,j), Y(i,j)];
        h_dir_pos = spatial_dir_channel_to_point(pos,scene,M,alpha_bt,dmin_bs);
        P_rad_NoRIS(i,j) = norm(h_dir_pos.'*wr_opt,2)^2;
        P_total_NoRIS(i,j) = abs(h_dir_pos.'*w_total_noris)^2;
    end
end

P_rad_NoRIS_dB = 10*log10(P_rad_NoRIS + 1e-13);
P_total_NoRIS_dB = 10*log10(P_total_NoRIS + 1e-13);

fig1 = figure('Color','w');
set(fig1,'Position',[120 120 680 560]);
imagesc(x_range,y_range,P_rad_NoRIS_dB);
set(gca,'YDir','normal');
axis equal tight;
colormap(gca,'jet');
colorbar;
hold on;
plot_scene_markers(scene);
title('Radar Power (No RIS)');
xlabel('x (m)'); ylabel('y (m)');
print(fig1, fullfile(outDir,'spatial_radar_power_noris.png'), '-dpng','-r300');
savefig(fig1, fullfile(outDir,'spatial_radar_power_noris.fig'));

fig2 = figure('Color','w');
set(fig2,'Position',[160 160 680 560]);
imagesc(x_range,y_range,P_total_NoRIS_dB);
set(gca,'YDir','normal');
axis equal tight;
colormap(gca,'jet');
colorbar;
hold on;
plot_scene_markers(scene);
title('ISAC Total Power (No RIS): wr + Wc');
xlabel('x (m)'); ylabel('y (m)');
print(fig2, fullfile(outDir,'spatial_total_power_noris.png'), '-dpng','-r300');
savefig(fig2, fullfile(outDir,'spatial_total_power_noris.fig'));

fprintf('[6/8] 生成No-RIS角域方向图...\n');
thetaScanDeg = -90:0.2:90;
thetaScanRad = thetaScanDeg*pi/180;
aScanBS = exp(-1j*(0:M-1)'*pi*sin(thetaScanRad))/sqrt(M);
ampBS = norm(Channel.hdt);
gain_noris = zeros(size(thetaScanDeg));
for t = 1:numel(thetaScanDeg)
    h_dir_t = ampBS*aScanBS(:,t);
    gain_noris(t) = norm(h_dir_t.'*wr_opt,2)^2;
end
gain_noris_dB = 10*log10(gain_noris + 1e-13);

fig3 = figure('Color','w');
plot(thetaScanDeg,gain_noris_dB,'-','LineWidth',1.6,'Color',[0.10 0.30 0.80]);
grid on;
xlabel('Angle (deg)');
ylabel('Radar beam gain (dB)');
title('No-RIS Radar Beampattern');
print(fig3, fullfile(outDir,'beampattern_noris.png'), '-dpng','-r300');
savefig(fig3, fullfile(outDir,'beampattern_noris.fig'));

fprintf('[7/8] 生成No-RIS eta权衡图...\n');
etas = 0:0.05:1;
SR_eta = zeros(size(etas));
SNR_eta = zeros(size(etas));
for i = 1:numel(etas)
    [Wc_i,wr_i] = design_w(Hk_noris,Channel.hdt,0,zeros(1,M),1,P,K,etas(i));
    SR_eta(i) = sum_rate(Hk_noris,Wc_i,sigmak2,wr_i,false);
    SNR_eta(i) = radar_snr_noris(Channel.hdt,wr_i,L,sigmat2,sigmar2,Wc_i,false);
end

fig4 = figure('Color','w');
set(fig4,'Position',[200 120 640 760]);
subplot(2,1,1);
plot(etas,10*log10(SNR_eta + 1e-13),'-','LineWidth',1.5,'Color',[0.80,0.10,0.10]);
hold on; grid on;
yline(10*log10(gammat),'--k','7 dB threshold','LabelHorizontalAlignment','left');
xlabel('\eta'); ylabel('Radar SNR (dB)');
title('No-RIS eta tradeoff');
subplot(2,1,2);
plot(etas,SR_eta,'-','LineWidth',1.5,'Color',[0.10,0.50,0.10]);
grid on;
xlabel('\eta'); ylabel('Sum-rate (bps/Hz)');
print(fig4, fullfile(outDir,'eta_tradeoff_noris.png'), '-dpng','-r300');
savefig(fig4, fullfile(outDir,'eta_tradeoff_noris.fig'));

fprintf('[8/8] 保存数据并输出汇总...\n');
save(fullfile(outDir,'noris_full_data.mat'), ...
    'M','K','L','P','gammat','x_range','y_range','X','Y', ...
    'P_rad_NoRIS','P_total_NoRIS','P_rad_NoRIS_dB','P_total_NoRIS_dB', ...
    'thetaScanDeg','gain_noris','gain_noris_dB', ...
    'etas','SR_eta','SNR_eta', ...
    'SR','gammat_opt','SNR_legacy','ptotal','Vsr_opt','scene');

fprintf('\n===== No-RIS Full Result Summary =====\n');
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
fprintf('Output directory      : %s\n', outDir);
fprintf('======================================\n\n');
end

function scene = build_paper_scene(K)
dg = 30;
scene.ris = [0, 0];
scene.bs = [-dg/sqrt(2), -dg/sqrt(2)];
scene.target = [3/sqrt(2), -3/sqrt(2)];

baseUsers = [
    -6, -2;
     6, -2
];

if K <= size(baseUsers,1)
    scene.users = baseUsers(1:K,:);
else
    extra = repmat(baseUsers(end,:), K-size(baseUsers,1), 1);
    scene.users = [baseUsers; extra];
end
end

function h_dir_pos = spatial_dir_channel_to_point(pos,scene,M,alpha_dir,dmin_bs)
vec_bs = pos - scene.bs;
d_bs = max(norm(vec_bs),dmin_bs);
theta_bs = atan2(vec_bs(2),vec_bs(1));

h_dir_pos = sqrt(1e-3*d_bs^(-alpha_dir)) * ...
    exp(-1j*(0:M-1)'*pi*sin(theta_bs))/sqrt(M);
end

function plot_scene_markers(scene)
plot(scene.bs(1),scene.bs(2),'ks','MarkerSize',7,'MarkerFaceColor','k');
plot(scene.ris(1),scene.ris(2),'kd','MarkerSize',7,'MarkerFaceColor',[0.2 0.2 0.2]);
plot(scene.target(1),scene.target(2),'rp','MarkerSize',12,'MarkerFaceColor','r');
plot(scene.users(:,1),scene.users(:,2),'bo','MarkerSize',6,'MarkerFaceColor','b');
legend('BS','RIS','Target','Users','Location','bestoutside');
end