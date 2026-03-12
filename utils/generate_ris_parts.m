function Channel = generate_ris_parts(N,baseline,varargin)
% generate_ris_parts - 生成RIS相关信道分量
%   Channel = generate_ris_parts(N,baseline) 按当前N随机生成RIS相关信道
%   Channel = generate_ris_parts(N,baseline,master) 从max-N主样本切片得到前N个单元
%
% 说明：
%   - 默认行为与历史版本一致：每次调用都会重新随机生成G/Hru
%   - 传入master后可实现“嵌套固定样本”比较：不同N共享同一底层随机实现

if ~isempty(varargin)
    master = varargin{1};
    assert(isfield(master,'hrt') && isfield(master,'G') && isfield(master,'Hru'), ...
        'master必须包含字段 hrt/G/Hru');
    assert(size(master.G,1) >= N && length(master.hrt) >= N && size(master.Hru,2) >= N, ...
        'master维度不足，无法切片到当前N');

    Channel.hrt = master.hrt(1:N,:);
    Channel.G = master.G(1:N,:);
    Channel.Hru = master.Hru(:,1:N);
    return;
end

M = baseline.M; K = baseline.K;
drt = baseline.drt; dg = baseline.dg; drk = baseline.drk;
alpha_rt = baseline.alpha_rt; alpha_g = baseline.alpha_g;
theta2 = baseline.theta2; thetar = baseline.thetar; kappa = baseline.kappa;
theta_ru = baseline.theta_ru;
hrt = sqrt(10^(-3)*drt^(-alpha_rt))*exp(-1j*(0:1:N-1)'*pi*sin(theta2));
GLos = sqrt(kappa/(1+kappa))*sqrt(10^(-3)*dg^(-alpha_g))*exp(-1j*(0:1:N-1)'*pi*sin(thetar))*exp(-1j*(0:1:M-1)*pi*sin(-thetar))/sqrt(M);
G = GLos + sqrt(1/(1+kappa))*sqrt(10^(-3)*dg^(-alpha_g))*(randn(N,M)+1i*randn(N,M))/sqrt(2*M);
Hru = zeros(K,N);
for k = 1:K
    theta_rk = theta_ru(k);
    Hru(k,:) = sqrt(10^(-3)*drk^(-baseline.alpha_rk))*(sqrt(kappa/(1+kappa))*exp(-1j*(0:1:N-1)*pi*sin(theta_rk)) + sqrt(1/(1+kappa))*(randn(1,N)+1i*randn(1,N))/sqrt(2));
end
Channel.hrt = hrt; Channel.G = G; Channel.Hru = Hru;
end