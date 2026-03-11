function baseline = generate_baseline(M,K)
drt = 3; dg = 50; drk = 8;
alpha_t = 2.4; alpha_k = 3.5;
theta2 = pi/4; thetar = pi/4;
theta1 = atan((dg*sin(thetar)-drt*cos(theta2))/(dg*cos(thetar)+drt*sin(theta2)));
dt = (dg*sin(thetar)-drt*cos(theta2))/sin(theta1);
hdt = sqrt(10^(-3)*dt^(-alpha_t))*exp(-1j*(0:1:M-1)'*pi*sin(theta1))/sqrt(M);
kappa = 10^(3/10);
theta_ru = pi/2*rand(K,1);
Hu = zeros(K,M);
for k = 1:K
    theta_rk = theta_ru(k);
    theta_k = atan((dg*sin(thetar)-drk*cos(theta_rk))/(dg*cos(thetar)+drk*sin(theta_rk)));
    dk = (dg*sin(thetar)-drk*cos(theta_rk))/sin(theta_k);
    Hu(k,:) = sqrt(10^(-3)*dk^(-alpha_k))*(sqrt(kappa/(1+kappa))*exp(-1j*(0:1:M-1)*pi*sin(theta_k))/sqrt(M) + sqrt(1/(1+kappa))*(randn(1,M)+1i*randn(1,M))/sqrt(2*M));
end
baseline.hdt = hdt;
baseline.theta_ru = theta_ru;
baseline.Hu = Hu;
baseline.M = M;
baseline.K = K;
baseline.drt = drt; baseline.dg = dg; baseline.drk = drk;
baseline.alpha_rt = 2.2; baseline.alpha_rk = 2.3; baseline.alpha_g = 2.2;
baseline.thetar = thetar; baseline.theta2 = theta2;
baseline.kappa = kappa;
end