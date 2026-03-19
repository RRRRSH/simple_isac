function [out, current_state] = evaluate_sr_baselines(Channel,baseline,P,K,L,sigmat2,sigmar2,sigmak2,gammat,varargin)
% evaluate_sr_baselines - 统一评估四种核心方案的通信总速率与雷达感知性能
% 支持温启动：通过传入 prev_state 加速收敛并保持单调性

prev_state = [];
if numel(varargin) >= 1 && (isstruct(varargin{1}) || isempty(varargin{1}))
    prev_state = varargin{1};
    varargin = varargin(2:end);
end

includeRadarInterferenceInRate = false;
includeCommInterferenceInRadar = false;
if numel(varargin) >= 1
    includeRadarInterferenceInRate = logical(varargin{1});
end
if numel(varargin) >= 2
    includeCommInterferenceInRadar = logical(varargin{2});
end

out = struct();
out.rate = struct('ris_joint',NaN,'ris_fixed',NaN,'ris_random',NaN,'no_ris',NaN);
out.snr = struct('ris_joint',NaN,'ris_fixed',NaN,'ris_random',NaN,'no_ris',NaN);
out.feasible = struct('ris_joint',false,'ris_fixed',false,'ris_random',false,'no_ris',false);

current_state = struct();

M = size(Channel.Hu,2);

% === 1. 先评估 RIS-Fixed (用作 Joint 的天然初值保障) ===
fprintf('  > 评估 RIS-Fixed (固定相位)...\n');
phi_fixed = compute_phi('fixed_target',Channel);
W0_fixed = [];
if ~isempty(prev_state) && isfield(prev_state,'W_fixed')
    W0_fixed = prev_state.W_fixed;
end
[Wc_fixed,wr_fixed,~,~,gammat_fixed_opt] = optimize_w_for_fixed_phi(Channel,phi_fixed,P,K,L,sigmat2,sigmar2,sigmak2,gammat,W0_fixed);
Hk_fixed = Channel.Hu + Channel.Hru*diag(phi_fixed)*Channel.G;
current_state.W_fixed = [Wc_fixed, wr_fixed];
out.rate.ris_fixed = sum_rate(Hk_fixed,Wc_fixed,sigmak2,wr_fixed,includeRadarInterferenceInRate);
out.snr.ris_fixed = gammat_fixed_opt;
out.feasible.ris_fixed = (out.snr.ris_fixed + 1e-12 >= gammat) && ~isnan(out.rate.ris_fixed);
if ~out.feasible.ris_fixed, out.rate.ris_fixed = NaN; end

% === 2. 评估 RIS-Joint (Proposed 采用 Warm-Start) ===
fprintf('  > 评估 RIS-Joint (Proposed)...\n');
if ~isempty(prev_state) && isfield(prev_state,'phi_joint') && isfield(prev_state,'W_joint')
    % 绝对温启动：继承上一个 Gamma_t 物理限制点留下的最优解
    phi0_joint = prev_state.phi_joint;
    W0_joint = prev_state.W_joint;
else
    % 冷启动：直接白嫖刚刚算完的 RIS-Fixed 输出充当优质物理基带，砍掉内置初值生成
    phi0_joint = phi_fixed;
    W0_joint = current_state.W_fixed; 
end
[phi_joint,Wc_joint,Wr_joint,~,gammat_joint] = run_joint_snr_optimization(Channel,P,K,L,sigmat2,sigmar2,sigmak2,gammat,phi0_joint,W0_joint);
Hk_joint = Channel.Hu + Channel.Hru*diag(phi_joint)*Channel.G;
rate_joint = sum_rate(Hk_joint,Wc_joint,sigmak2,Wr_joint(:,K+1:end),includeRadarInterferenceInRate);
feasible_joint = (gammat_joint + 1e-12 >= gammat) && ~isnan(rate_joint);

% === 2b. 自愈机制：如果温启动导致无解，立刻尝试用刚算的 Fixed 结果做冷启动 ===
if ~feasible_joint
    fprintf('        -> [自愈] 温启动失败(可能门限过高)，尝试使用 Fixed 结果进行冷启动修复...\n');
    phi0_cold = phi_fixed;
    W0_cold = current_state.W_fixed;
    [phi_joint,Wc_joint,Wr_joint,~,gammat_joint] = run_joint_snr_optimization(Channel,P,K,L,sigmat2,sigmar2,sigmak2,gammat,phi0_cold,W0_cold);
    Hk_joint = Channel.Hu + Channel.Hru*diag(phi_joint)*Channel.G;
    rate_joint = sum_rate(Hk_joint,Wc_joint,sigmak2,Wr_joint(:,K+1:end),includeRadarInterferenceInRate);
    feasible_joint = (gammat_joint + 1e-12 >= gammat) && ~isnan(rate_joint);
end

current_state.phi_joint = phi_joint;
current_state.W_joint = Wr_joint;
out.rate.ris_joint = rate_joint;
out.snr.ris_joint = gammat_joint;
out.feasible.ris_joint = feasible_joint;
if ~out.feasible.ris_joint, out.rate.ris_joint = NaN; end

% === 3. 评估 RIS-Random ===
fprintf('  > 评估 RIS-Random (随机相位)...\n');
% 维持随机初值的生命周期，确保跨步长扫描时是同一套障碍空间环境配置
if ~isempty(prev_state) && isfield(prev_state,'phi_random')
    phi_random = prev_state.phi_random;
else
    phi_random = compute_phi('random',Channel);
end
W0_random = [];
if ~isempty(prev_state) && isfield(prev_state,'W_random')
    W0_random = prev_state.W_random;
end
[Wc_random,wr_random,~,~,gammat_random_opt] = optimize_w_for_fixed_phi(Channel,phi_random,P,K,L,sigmat2,sigmar2,sigmak2,gammat,W0_random);
Hk_random = Channel.Hu + Channel.Hru*diag(phi_random)*Channel.G;
current_state.phi_random = phi_random;
current_state.W_random = [Wc_random, wr_random];
out.rate.ris_random = sum_rate(Hk_random,Wc_random,sigmak2,wr_random,includeRadarInterferenceInRate);
out.snr.ris_random = gammat_random_opt;
out.feasible.ris_random = (out.snr.ris_random + 1e-12 >= gammat) && ~isnan(out.rate.ris_random);
if ~out.feasible.ris_random, out.rate.ris_random = NaN; end

% === 4. 评估 No-RIS ===
fprintf('  > 评估 No-RIS (纯基站)...\n');
Channel_nr.hdt = baseline.hdt; Channel_nr.Hu = baseline.Hu;
Channel_nr.hrt = 0; Channel_nr.G = zeros(1,M); Channel_nr.Hru = zeros(K,1);
phi_nr = 1;
W0_nr = [];
if ~isempty(prev_state) && isfield(prev_state,'W_noris')
    W0_nr = prev_state.W_noris;
end
[Wc_nr,wr_nr,~,~,gammat_nr_opt] = optimize_w_for_fixed_phi(Channel_nr,phi_nr,P,K,L,sigmat2,sigmar2,sigmak2,gammat,W0_nr);
Hk_nr = baseline.Hu;
current_state.W_noris = [Wc_nr, wr_nr];
out.rate.no_ris = sum_rate(Hk_nr,Wc_nr,sigmak2,wr_nr,includeRadarInterferenceInRate);
out.snr.no_ris = gammat_nr_opt;
out.feasible.no_ris = (out.snr.no_ris + 1e-12 >= gammat) && ~isnan(out.rate.no_ris);
if ~out.feasible.no_ris, out.rate.no_ris = NaN; end

end