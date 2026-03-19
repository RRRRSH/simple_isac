function plot_spatial_heatmap_nsp_vs_mrt(outDir,gammatFixed_dB)
% plot_spatial_heatmap_nsp_vs_mrt
% 绘制固定相位下 No-RIS / RIS+MRT / RIS+NSP 的空间热力图。

baseDir = fileparts(mfilename('fullpath'));
addpath(fullfile(baseDir,'utils'));
rng(1);

if nargin < 1 || isempty(outDir)
    outRoot = fullfile(baseDir,'out');
    if ~exist(outRoot,'dir'), mkdir(outRoot); end
    runStamp = datestr(now,'yyyymmdd_HHMMSS');
    outDir = fullfile(outRoot,['plot_spatial_heatmap_nsp_vs_mrt_' runStamp]);
end
if nargin < 2 || isempty(gammatFixed_dB)
    gammatFixed_dB = 7;
end
if ~exist(outDir,'dir'), mkdir(outDir); end
fprintf('Heatmap outputs will be saved to: %s\n', outDir);

% System settings
M = 6;
K = 4;
N = 64;
L = 1024;
P = 10.^((30-30)/10);  % 30 dBm -> 1W
gammat = 10^(gammatFixed_dB/10);
sigmar2 = 1e-12;
sigmat2 = 1;
includeCommInterferenceInRadar = false;

% Geometry and channel
scene = build_paper_scene(K);
baseline = generate_baseline(M,K);
ris = generate_ris_parts(N,baseline);

Channel.hdt = baseline.hdt;
Channel.Hu = baseline.Hu;
Channel.hrt = ris.hrt;
Channel.G = ris.G;
Channel.Hru = ris.Hru;
phi_fixed = compute_phi('fixed_target',Channel);
Hk_ris = Channel.Hu + Channel.Hru*diag(phi_fixed)*Channel.G;
sigmak2 = 1e-12;

% RIS + 优化波束（默认NSP初值）
[Wc_nsp,wr_nsp] = optimize_w_for_fixed_phi(Channel,phi_fixed,P,K,L,sigmat2,sigmar2,sigmak2,gammat);

% RIS + 优化波束（MRT初值）
eta0_mrt = 0.5;
[Wc0_mrt,wr0_mrt] = design_w(Hk_ris,Channel.hdt,Channel.hrt,Channel.G,phi_fixed,P,K,eta0_mrt,'mrt');
W0_mrt = [Wc0_mrt, wr0_mrt, zeros(M,max(M-1,0))];
Prms.M = M; Prms.N = N; Prms.K = K;
Prms.sigmar2 = sigmar2; Prms.sigmak2 = sigmak2; Prms.sigmat2 = sigmat2;
Prms.Nmax = 40; Prms.res_th = 1e-3; Prms.L = L;
Prms.gammat = gammat; Prms.P = P; Prms.freezePhi = true;
try
    [W_mrt,~,~,~] = get_W_phi_SNR(Prms,Channel,phi_fixed,W0_mrt);
catch
    W_mrt = W0_mrt;
end
Wc_mrt = W_mrt(:,1:K);
wr_mrt = W_mrt(:,K+1:end);

% No-RIS baseline: 通过退化信道使用优化器
Channel_nr.hdt = baseline.hdt;
Channel_nr.Hu = baseline.Hu;
Channel_nr.hrt = 0;
Channel_nr.G = zeros(1,M);
Channel_nr.Hru = zeros(K,1);
phi_nr = 1;
[Wc_nr,wr_nr] = optimize_w_for_fixed_phi(Channel_nr,phi_nr,P,K,L,sigmat2,sigmar2,sigmak2,gammat);

allPts = [scene.bs; scene.ris; scene.target; scene.users];
margin = 10;
x_range = floor(min(allPts(:,1))-margin):0.5:ceil(max(allPts(:,1))+margin);
y_range = floor(min(allPts(:,2))-margin):0.5:ceil(max(allPts(:,2))+margin);
[X,Y] = meshgrid(x_range,y_range);

alpha_bt = 2.2;
alpha_rt = 2.2;
dmin_bs = 1e-6;
dmin_ris = 1e-6;

P_rad_nr = zeros(size(X));
P_rad_mrt = zeros(size(X));
P_rad_nsp = zeros(size(X));

w_total_nr = wr_nr + sum(Wc_nr,2);
w_total_mrt = wr_mrt + sum(Wc_mrt,2);
w_total_nsp = wr_nsp + sum(Wc_nsp,2);

P_total_nr = zeros(size(X));
P_total_mrt = zeros(size(X));
P_total_nsp = zeros(size(X));

for i = 1:size(X,1)
    for j = 1:size(X,2)
        pos = [X(i,j), Y(i,j)];
        [h_dir_pos, h_ris_pos] = spatial_channels_to_point(pos,scene,M,N,alpha_bt,alpha_rt,dmin_bs,dmin_ris);
        h_eff_nr = h_dir_pos.';
        h_eff_ris = h_dir_pos.' + h_ris_pos.'*diag(phi_fixed)*Channel.G;

        P_rad_nr(i,j) = abs(h_eff_nr*wr_nr)^2;
        P_rad_mrt(i,j) = abs(h_eff_ris*wr_mrt)^2;
        P_rad_nsp(i,j) = abs(h_eff_ris*wr_nsp)^2;

        P_total_nr(i,j) = abs(h_eff_nr*w_total_nr)^2;
        P_total_mrt(i,j) = abs(h_eff_ris*w_total_mrt)^2;
        P_total_nsp(i,j) = abs(h_eff_ris*w_total_nsp)^2;
    end
end

P_rad_nr_dB = 10*log10(P_rad_nr + 1e-13);
P_rad_mrt_dB = 10*log10(P_rad_mrt + 1e-13);
P_rad_nsp_dB = 10*log10(P_rad_nsp + 1e-13);
P_total_nr_dB = 10*log10(P_total_nr + 1e-13);
P_total_mrt_dB = 10*log10(P_total_mrt + 1e-13);
P_total_nsp_dB = 10*log10(P_total_nsp + 1e-13);

cLimRadar = [min([P_rad_nr_dB(:);P_rad_mrt_dB(:);P_rad_nsp_dB(:)]), max([P_rad_nr_dB(:);P_rad_mrt_dB(:);P_rad_nsp_dB(:)])];
cLimTotal = [min([P_total_nr_dB(:);P_total_mrt_dB(:);P_total_nsp_dB(:)]), max([P_total_nr_dB(:);P_total_mrt_dB(:);P_total_nsp_dB(:)])];

fig = figure('Color','w');
set(fig,'Position',[80 80 1500 860]);

subplot(2,3,1); imagesc(x_range,y_range,P_rad_nr_dB); set(gca,'YDir','normal'); axis equal tight; caxis(cLimRadar); colormap(gca,'jet'); colorbar; hold on; plot_scene_markers(scene); title('Radar Power: No-RIS'); xlabel('x (m)'); ylabel('y (m)');
subplot(2,3,2); imagesc(x_range,y_range,P_rad_mrt_dB); set(gca,'YDir','normal'); axis equal tight; caxis(cLimRadar); colormap(gca,'jet'); colorbar; hold on; plot_scene_markers(scene); title('Radar Power: RIS+MRT'); xlabel('x (m)'); ylabel('y (m)');
subplot(2,3,3); imagesc(x_range,y_range,P_rad_nsp_dB); set(gca,'YDir','normal'); axis equal tight; caxis(cLimRadar); colormap(gca,'jet'); colorbar; hold on; plot_scene_markers(scene); title('Radar Power: RIS+NSP'); xlabel('x (m)'); ylabel('y (m)');

subplot(2,3,4); imagesc(x_range,y_range,P_total_nr_dB); set(gca,'YDir','normal'); axis equal tight; caxis(cLimTotal); colormap(gca,'jet'); colorbar; hold on; plot_scene_markers(scene); title('Total Power: No-RIS'); xlabel('x (m)'); ylabel('y (m)');
subplot(2,3,5); imagesc(x_range,y_range,P_total_mrt_dB); set(gca,'YDir','normal'); axis equal tight; caxis(cLimTotal); colormap(gca,'jet'); colorbar; hold on; plot_scene_markers(scene); title('Total Power: RIS+MRT'); xlabel('x (m)'); ylabel('y (m)');
subplot(2,3,6); imagesc(x_range,y_range,P_total_nsp_dB); set(gca,'YDir','normal'); axis equal tight; caxis(cLimTotal); colormap(gca,'jet'); colorbar; hold on; plot_scene_markers(scene); title('Total Power: RIS+NSP'); xlabel('x (m)'); ylabel('y (m)');

print(fig, fullfile(outDir,'spatial_heatmap_noris_mrt_nsp.png'), '-dpng','-r300');
savefig(fig, fullfile(outDir,'spatial_heatmap_noris_mrt_nsp.fig'));

save(fullfile(outDir,'spatial_heatmap_noris_mrt_nsp_data.mat'), ...
    'x_range','y_range','X','Y', ...
    'P_rad_nr','P_rad_mrt','P_rad_nsp','P_total_nr','P_total_mrt','P_total_nsp', ...
    'P_rad_nr_dB','P_rad_mrt_dB','P_rad_nsp_dB','P_total_nr_dB','P_total_mrt_dB','P_total_nsp_dB', ...
    'eta_nr','eta_mrt','eta_nsp','gammat','gammatFixed_dB','scene');
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

function [h_dir_pos, h_ris_pos] = spatial_channels_to_point(pos,scene,M,N,alpha_dir,alpha_ris,dmin_bs,dmin_ris)
vec_bs = pos - scene.bs;
d_bs = max(norm(vec_bs),dmin_bs);
theta_bs = atan2(vec_bs(2),vec_bs(1));

vec_ris = pos - scene.ris;
d_ris = max(norm(vec_ris),dmin_ris);
theta_ris = atan2(vec_ris(2),vec_ris(1));

h_dir_pos = sqrt(1e-3*d_bs^(-alpha_dir)) * ...
    exp(-1j*(0:M-1)'*pi*sin(theta_bs))/sqrt(M);

h_ris_pos = sqrt(1e-3*d_ris^(-alpha_ris)) * ...
    exp(-1j*(0:N-1)'*pi*sin(theta_ris));
end

function plot_scene_markers(scene)
plot(scene.bs(1),scene.bs(2),'ks','MarkerSize',6,'MarkerFaceColor','k');
plot(scene.ris(1),scene.ris(2),'kd','MarkerSize',6,'MarkerFaceColor',[0.2 0.2 0.2]);
plot(scene.target(1),scene.target(2),'rp','MarkerSize',10,'MarkerFaceColor','r');
plot(scene.users(:,1),scene.users(:,2),'bo','MarkerSize',5,'MarkerFaceColor','b');
end
