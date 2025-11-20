function Channel = generate_ris_parts(N,baseline)
M = baseline.M; K = baseline.K;
drt = baseline.drt; dg = baseline.dg; drk = baseline.drk;
alpha_rt = baseline.alpha_rt; alpha_rk = baseline.alpha_rk; alpha_g = baseline.alpha_g;
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