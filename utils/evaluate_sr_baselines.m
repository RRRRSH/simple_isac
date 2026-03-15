function out = evaluate_sr_baselines(Channel,baseline,P,K,L,sigmat2,sigmar2,sigmak2,gammat,varargin)
% evaluate_sr_baselines - 统一评估四种方案的通信总速率
% 方案: Proposed(联合优化), Separate(固定相位+波束优化), Comm-only(纯通信), No-RIS

includeRadarInterferenceInRate = false;
includeCommInterferenceInRadar = false;
runMode = 'all';
if numel(varargin) >= 1
    includeRadarInterferenceInRate = logical(varargin{1});
end
if numel(varargin) >= 2
    includeCommInterferenceInRadar = logical(varargin{2});
end
if numel(varargin) >= 3
    runMode = lower(string(varargin{3}));
end

isNoRISOnly = any(runMode == ["noris_only","no_ris_only"]);
isNoRISVsFixed = any(runMode == ["noris_vs_fixed","no_ris_vs_fixed"]);
isNoRISVsFixedNSPMRT = any(runMode == ["noris_vs_fixed_nsp_vs_mrt","no_ris_vs_fixed_nsp_vs_mrt"]);

out = struct();
out.rate = struct('proposed',NaN,'separate',NaN,'comm_only',NaN,'no_ris',NaN);
out.snr = struct('proposed',NaN,'separate',NaN,'comm_only',NaN,'no_ris',NaN);
out.feasible = struct('proposed',false,'separate',false,'comm_only',true,'no_ris',false);
out.rate.fixed_nsp = NaN;
out.rate.fixed_mrt = NaN;
out.snr.fixed_nsp = NaN;
out.snr.fixed_mrt = NaN;
out.feasible.fixed_nsp = false;
out.feasible.fixed_mrt = false;
out.meta = struct('runMode',char(runMode));

if ~(isNoRISOnly || isNoRISVsFixed || isNoRISVsFixedNSPMRT)
    % Proposed: 联合优化（相位+波束）
    [phi_joint,Wc_joint,Wr_joint] = run_joint_snr_optimization(Channel,P,K,L,sigmat2,sigmar2,sigmak2,gammat);
    Hk_joint = Channel.Hu + Channel.Hru*diag(phi_joint)*Channel.G;
    Wr_radar_only = Wr_joint(:,K+1:end);
    out.rate.proposed = sum_rate(Hk_joint,Wc_joint,sigmak2,Wr_radar_only,includeRadarInterferenceInRate);
    out.snr.proposed = radar_snr(Channel.hdt,Channel.hrt,Channel.G,phi_joint,Wr_joint,L,sigmat2,sigmar2,Wc_joint,includeCommInterferenceInRadar);
    out.feasible.proposed = (out.snr.proposed + 1e-12 >= gammat);
    if ~out.feasible.proposed
        out.rate.proposed = NaN;
    end

    % Separate: 固定相位后仅优化波束
    phi_sep = compute_phi('fixed_target',Channel);
    [Wc_sep,wr_sep] = optimize_w_for_fixed_phi(Channel,phi_sep,P,K,L,sigmat2,sigmar2,sigmak2,gammat);
    Hk_sep = Channel.Hu + Channel.Hru*diag(phi_sep)*Channel.G;
    out.rate.separate = sum_rate(Hk_sep,Wc_sep,sigmak2,wr_sep,includeRadarInterferenceInRate);
    out.snr.separate = radar_snr(Channel.hdt,Channel.hrt,Channel.G,phi_sep,wr_sep,L,sigmat2,sigmar2,Wc_sep,includeCommInterferenceInRadar);
    out.feasible.separate = (out.snr.separate + 1e-12 >= gammat);
    if ~out.feasible.separate
        out.rate.separate = NaN;
    end

    % Comm-only: 不加感知约束，eta=0，按sum-rate目标扫描相位
    phi_comm = scan_phi('sumrate',Channel,P,K,L,sigmat2,sigmar2,0,sigmak2, ...
        includeRadarInterferenceInRate,includeCommInterferenceInRadar);
    Hk_comm = Channel.Hu + Channel.Hru*diag(phi_comm)*Channel.G;
    [Wc_comm,~] = design_w(Hk_comm,Channel.hdt,Channel.hrt,Channel.G,phi_comm,P,K,0);
    out.rate.comm_only = sum_rate(Hk_comm,Wc_comm,sigmak2);
    out.snr.comm_only = radar_snr(Channel.hdt,Channel.hrt,Channel.G,phi_comm,zeros(size(Channel.hdt)), ...
        L,sigmat2,sigmar2,Wc_comm,includeCommInterferenceInRadar);
elseif isNoRISOnly
    % 快速模式：其他三方案不计算，按0填充便于同一数据结构出图
    out.rate.proposed = 0;
    out.rate.separate = 0;
    out.rate.comm_only = 0;
    out.snr.proposed = 0;
    out.snr.separate = 0;
    out.snr.comm_only = 0;
elseif isNoRISVsFixed
    % 双线模式：仅计算固定相位与No-RIS，跳过其余两方案
    out.rate.proposed = 0;
    out.rate.comm_only = 0;
    out.snr.proposed = 0;
    out.snr.comm_only = 0;

    phi_sep = compute_phi('fixed_target',Channel);
    % 为保证与可行性判据一致，双线快速模式采用同口径启发式固定相位波束
    eta_sep = select_eta(gammat,Channel,phi_sep,P,K,L,sigmat2,sigmar2,includeCommInterferenceInRadar,'mrt');
    Hk_sep = Channel.Hu + Channel.Hru*diag(phi_sep)*Channel.G;
    [Wc_sep,wr_sep] = design_w(Hk_sep,Channel.hdt,Channel.hrt,Channel.G,phi_sep,P,K,eta_sep,'mrt');
    out.rate.separate = sum_rate(Hk_sep,Wc_sep,sigmak2,wr_sep,includeRadarInterferenceInRate);
    out.snr.separate = radar_snr(Channel.hdt,Channel.hrt,Channel.G,phi_sep,wr_sep,L,sigmat2,sigmar2,Wc_sep,includeCommInterferenceInRadar);
    out.feasible.separate = (out.snr.separate + 1e-12 >= gammat);
    if ~out.feasible.separate
        out.rate.separate = NaN;
    end
elseif isNoRISVsFixedNSPMRT
    % 三线模式：No-RIS / RIS+MRT / RIS+NSP
    out.rate.proposed = 0;
    out.rate.comm_only = 0;
    out.rate.separate = 0;
    out.snr.proposed = 0;
    out.snr.comm_only = 0;
    out.snr.separate = 0;

    phi_sep = compute_phi('fixed_target',Channel);
    Hk_sep = Channel.Hu + Channel.Hru*diag(phi_sep)*Channel.G;

    eta_nsp = select_eta(gammat,Channel,phi_sep,P,K,L,sigmat2,sigmar2,includeCommInterferenceInRadar,'nsp');
    [Wc_nsp,wr_nsp] = design_w(Hk_sep,Channel.hdt,Channel.hrt,Channel.G,phi_sep,P,K,eta_nsp,'nsp');
    out.rate.fixed_nsp = sum_rate(Hk_sep,Wc_nsp,sigmak2,wr_nsp,includeRadarInterferenceInRate);
    out.snr.fixed_nsp = radar_snr(Channel.hdt,Channel.hrt,Channel.G,phi_sep,wr_nsp,L,sigmat2,sigmar2,Wc_nsp,includeCommInterferenceInRadar);
    out.feasible.fixed_nsp = (out.snr.fixed_nsp + 1e-12 >= gammat);
    if ~out.feasible.fixed_nsp
        out.rate.fixed_nsp = NaN;
    end

    eta_mrt = select_eta(gammat,Channel,phi_sep,P,K,L,sigmat2,sigmar2,includeCommInterferenceInRadar,'mrt');
    [Wc_mrt,wr_mrt] = design_w(Hk_sep,Channel.hdt,Channel.hrt,Channel.G,phi_sep,P,K,eta_mrt,'mrt');
    out.rate.fixed_mrt = sum_rate(Hk_sep,Wc_mrt,sigmak2,wr_mrt,includeRadarInterferenceInRate);
    out.snr.fixed_mrt = radar_snr(Channel.hdt,Channel.hrt,Channel.G,phi_sep,wr_mrt,L,sigmat2,sigmar2,Wc_mrt,includeCommInterferenceInRadar);
    out.feasible.fixed_mrt = (out.snr.fixed_mrt + 1e-12 >= gammat);
    if ~out.feasible.fixed_mrt
        out.rate.fixed_mrt = NaN;
    end
end

% No-RIS: 仅基站波束，不使用RIS
Hk_nr = baseline.Hu;
eta_nr = 1;
for e = 0:0.05:1
    [Wc_tmp,wr_tmp] = design_w(Hk_nr,baseline.hdt,0,zeros(1,size(Hk_nr,2)),1,P,K,e,'mrt');
    snr_tmp = radar_snr_noris(baseline.hdt,wr_tmp,L,sigmat2,sigmar2,Wc_tmp,includeCommInterferenceInRadar);
    if snr_tmp >= gammat
        eta_nr = e;
        break;
    end
end
[Wc_nr,wr_nr] = design_w(Hk_nr,baseline.hdt,0,zeros(1,size(Hk_nr,2)),1,P,K,eta_nr,'mrt');
out.rate.no_ris = sum_rate(Hk_nr,Wc_nr,sigmak2,wr_nr,includeRadarInterferenceInRate);
out.snr.no_ris = radar_snr_noris(baseline.hdt,wr_nr,L,sigmat2,sigmar2,Wc_nr,includeCommInterferenceInRadar);
out.feasible.no_ris = (out.snr.no_ris + 1e-12 >= gammat);
if ~out.feasible.no_ris
    out.rate.no_ris = NaN;
end
end