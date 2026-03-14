function plot_spatial_heatmap()
% plot_spatial_heatmap
% Generate three core spatial heatmaps for RIS-ISAC visual evidence:
% 1) Radar absolute power comparison (No RIS vs RIS-Fixed)
% 2) Spatial power difference map Delta P (RIS - No RIS, in dB)
% 3) ISAC total power map under RIS-Fixed (radar + communication beams)

baseDir = fileparts(mfilename('fullpath'));
addpath(fullfile(baseDir,'utils'));
rng(1);

% Output directory with timestamp
outRoot = fullfile(baseDir,'out');
if ~exist(outRoot,'dir'), mkdir(outRoot); end
runStamp = datestr(now,'yyyymmdd_HHMMSS');
outDir = fullfile(outRoot,['plot_spatial_heatmap_' runStamp]);
if ~exist(outDir,'dir'), mkdir(outDir); end
fprintf('Results will be saved to: %s\n', outDir);

% System parameters (consistent with run_simple_isac)
M = 6;
K = 2;
L = 1024;
sigmar2 = 1e-12;
sigmat2 = 1;
P = 10.^(32/10-3);
N0 = 64;
alpha_bt = 2.2;   % BS->target / BS->RIS-related pathloss exponent
alpha_bu = 3.5;   % BS->user pathloss exponent
alpha_rt = 2.2;
alpha_rk = 2.3;
alpha_g = 2.2;
kappa = 10^(3/10);
includeCommInterferenceInRadar = false;
dmin_bs = 1e-6;
dmin_ris = 1e-6;

% Fixed color axes to avoid visual exaggeration across runs
cLimAbs = [-120, -40];
cLimDelta = [-10, 10];
cLimTotal = [-120, -30];

% Build paper-style geometry and generate channels consistently from it
scene = build_paper_scene(K);

[hdt,hrt] = spatial_channels_to_point(scene.target,scene,M,N0,alpha_bt,alpha_rt,dmin_bs,dmin_ris);
Hu = zeros(K,M);
Hru = zeros(K,N0);
for k = 1:K
    [hdu,hru] = spatial_channels_to_point(scene.users(k,:),scene,M,N0,alpha_bu,alpha_rk,dmin_bs,dmin_ris);
    Hu(k,:) = hdu.';
    Hru(k,:) = hru.';
end
G = bs_ris_channel(scene,M,N0,alpha_g,kappa);

Channel0.hdt = hdt;
Channel0.Hu = Hu;
Channel0.hrt = hrt;
Channel0.G = G;
Channel0.Hru = Hru;

% Fixed-target RIS phase and beams
phi_fixed = compute_phi('fixed_target',Channel0);
Hk_ris = Channel0.Hu + Channel0.Hru*diag(phi_fixed)*Channel0.G;

gammat_legacy = 10^0.7;
gammat = gammat_legacy*(norm(Channel0.hdt)^2);
eta_ris = select_eta(gammat,Channel0,phi_fixed,P,K,L,sigmat2,sigmar2,includeCommInterferenceInRadar);
[Wc_ris,wr_ris] = design_w(Hk_ris,Channel0.hdt,Channel0.hrt,Channel0.G,phi_fixed,P,K,eta_ris);

% No-RIS beams
Hk_noris = Channel0.Hu;
G_noris = zeros(1,M);
phi_noris = 1;
eta_noris = 0;
for e = 0:0.05:1
    [Wc_tmp,wr_tmp] = design_w(Hk_noris,Channel0.hdt,0,G_noris,phi_noris,P,K,e);
    if radar_snr_noris(Channel0.hdt,wr_tmp,L,sigmat2,sigmar2,Wc_tmp,includeCommInterferenceInRadar) >= gammat
        eta_noris = e;
        break;
    end
end
[~,wr_noris] = design_w(Hk_noris,Channel0.hdt,0,G_noris,phi_noris,P,K,eta_noris);

% 2D virtual space grid (fixed for reproducible visual comparison)
allPts = [scene.bs; scene.ris; scene.target; scene.users];
margin = 10;
x_range = floor(min(allPts(:,1))-margin):0.5:ceil(max(allPts(:,1))+margin);
y_range = floor(min(allPts(:,2))-margin):0.5:ceil(max(allPts(:,2))+margin);
[X,Y] = meshgrid(x_range,y_range);

P_rad_NoRIS = zeros(size(X));
P_rad_RIS = zeros(size(X));
P_total_RIS = zeros(size(X));

w_total_ris = wr_ris + sum(Wc_ris,2);

% Scan each grid point
for i = 1:size(X,1)
    for j = 1:size(X,2)
        pos = [X(i,j), Y(i,j)];
        [h_dir_pos, h_ris_pos] = spatial_channels_to_point(pos,scene,M,N0,alpha_bt,alpha_rt,dmin_bs,dmin_ris);

        % Effective channels
        h_eff_NoRIS = h_dir_pos.';
        h_eff_RIS = h_dir_pos.' + h_ris_pos.'*diag(phi_fixed)*Channel0.G;

        % Forward transmit spatial power (radiation map)
        P_rad_NoRIS(i,j) = abs(h_eff_NoRIS*wr_noris)^2;
        P_rad_RIS(i,j) = abs(h_eff_RIS*wr_ris)^2;
        P_total_RIS(i,j) = abs(h_eff_RIS*w_total_ris)^2;
    end
end

% dB conversion with floor buffer to avoid -Inf
P_rad_NoRIS_dB = 10*log10(P_rad_NoRIS + 1e-13);
P_rad_RIS_dB = 10*log10(P_rad_RIS + 1e-13);
P_total_RIS_dB = 10*log10(P_total_RIS + 1e-13);

% Difference map
Delta_P_dB = P_rad_RIS_dB - P_rad_NoRIS_dB;

% Figure 1: Radar absolute power comparison
fig1 = figure('Color','w');
set(fig1,'Position',[80 80 1200 500]);

subplot(1,2,1);
imagesc(x_range,y_range,P_rad_NoRIS_dB);
set(gca,'YDir','normal');
axis equal tight;
caxis(cLimAbs);
colormap(gca,'jet');
colorbar;
hold on;
plot_scene_markers(scene);
title('Radar Power (No RIS)');
xlabel('x (m)'); ylabel('y (m)');

subplot(1,2,2);
imagesc(x_range,y_range,P_rad_RIS_dB);
set(gca,'YDir','normal');
axis equal tight;
caxis(cLimAbs);
colormap(gca,'jet');
colorbar;
hold on;
plot_scene_markers(scene);
title('Radar Power (RIS-Fixed)');
xlabel('x (m)'); ylabel('y (m)');

print(fig1, fullfile(outDir,'spatial_radar_power_compare.png'), '-dpng','-r300');
savefig(fig1, fullfile(outDir,'spatial_radar_power_compare.fig'));

% Figure 2: Difference map Delta P = RIS - No RIS
fig2 = figure('Color','w');
set(fig2,'Position',[120 120 680 560]);
imagesc(x_range,y_range,Delta_P_dB);
set(gca,'YDir','normal');
axis equal tight;
caxis(cLimDelta);
colormap(gca,'jet');
cb2 = colorbar;
cb2.Label.String = '\DeltaP (dB)';
hold on;
plot_scene_markers(scene);
title('Spatial Power Gain: RIS-Fixed - No RIS');
xlabel('x (m)'); ylabel('y (m)');

print(fig2, fullfile(outDir,'spatial_delta_power_ris_minus_noris.png'), '-dpng','-r300');
savefig(fig2, fullfile(outDir,'spatial_delta_power_ris_minus_noris.fig'));

% Figure 3: ISAC total power map under RIS-Fixed
fig3 = figure('Color','w');
set(fig3,'Position',[160 160 680 560]);
imagesc(x_range,y_range,P_total_RIS_dB);
set(gca,'YDir','normal');
axis equal tight;
caxis(cLimTotal);
colormap(gca,'jet');
colorbar;
hold on;
plot_scene_markers(scene);
title('ISAC Total Power (RIS-Fixed): wr + Wc');
xlabel('x (m)'); ylabel('y (m)');

print(fig3, fullfile(outDir,'spatial_total_power_ris_fixed.png'), '-dpng','-r300');
savefig(fig3, fullfile(outDir,'spatial_total_power_ris_fixed.fig'));

% Save matrices for paper post-processing
save(fullfile(outDir,'spatial_heatmap_data.mat'), ...
    'x_range','y_range','X','Y', ...
    'P_rad_NoRIS','P_rad_RIS','P_total_RIS', ...
    'P_rad_NoRIS_dB','P_rad_RIS_dB','P_total_RIS_dB','Delta_P_dB', ...
    'scene','eta_noris','eta_ris','cLimAbs','cLimDelta','cLimTotal');
end

function scene = build_paper_scene(K)
% User-specified clean layout for spatial visualization
% RIS at origin, BS at [-dg/sqrt(2), -dg/sqrt(2)] with dg=30,
% two users at [-6,-2], [6,-2], target at distance 3 and angle -45 deg.

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
% Generate BS->point and RIS->point channels for each spatial sample point.

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

function G = bs_ris_channel(scene,M,N,alpha_g,kappa)
% RIS-BS channel with LOS + weak scattering to avoid rank deficiency.
vec = scene.ris - scene.bs;
dg = max(norm(vec),1e-6);
theta = atan2(vec(2),vec(1));

GLos = sqrt(kappa/(1+kappa))*sqrt(1e-3*dg^(-alpha_g)) * ...
    exp(-1j*(0:N-1)'*pi*sin(theta)) * exp(-1j*(0:M-1)*pi*sin(-theta))/sqrt(M);
GNlos = sqrt(1/(1+kappa))*sqrt(1e-3*dg^(-alpha_g))*(randn(N,M)+1i*randn(N,M))/sqrt(2*M);
G = GLos + GNlos;
end

function plot_scene_markers(scene)
plot(scene.bs(1),scene.bs(2),'ks','MarkerSize',7,'MarkerFaceColor','k');
plot(scene.ris(1),scene.ris(2),'kd','MarkerSize',7,'MarkerFaceColor',[0.2 0.2 0.2]);
plot(scene.target(1),scene.target(2),'rp','MarkerSize',12,'MarkerFaceColor','r');
plot(scene.users(:,1),scene.users(:,2),'bo','MarkerSize',6,'MarkerFaceColor','b');
legend('BS','RIS','Target','Users','Location','bestoutside');
end
