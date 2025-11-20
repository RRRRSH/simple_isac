% 函数：generate_channels
% 功能：为RIS辅助的ISAC系统生成所有必要的信道矩阵
% 输入参数：
%   M - 基站天线数量
%   N - RIS单元数量
%   K - 通信用户数量
% 输出参数：
%   Channel - 包含所有信道矩阵的结构体
function Channel = generate_channels(M,N,K)
    % 系统几何参数设置
    drt = 3;          % 目标到RIS的距离(单位：米)
    dg = 50;          % RIS到基站的距离(单位：米)
    drk = 8;          % 用户到RIS的距离(单位：米)
    
    % 路径损耗指数设置
    alpha_t = 2.4;    % 直接雷达信道(目标-基站)的路径损耗指数
    alpha_rt = 2.2;   % RIS-目标信道的路径损耗指数
    alpha_k = 3.5;    % 直接用户信道(用户-基站)的路径损耗指数
    alpha_rk = 2.3;   % RIS-用户信道的路径损耗指数
    alpha_g = 2.2;    % RIS-基站信道的路径损耗指数
    
    % 角度参数设置
    theta2 = pi/4;    % 目标到RIS的入射角(弧度)
    thetar = pi/4;    % RIS到基站的入射角(弧度)
    
    % 计算目标到基站的入射角(基于几何关系)
    theta1 = atan((dg*sin(thetar)-drt*cos(theta2))/(dg*cos(thetar)+drt*sin(theta2)));
    
    % 计算目标到基站的直接距离
    dt = (dg*sin(thetar)-drt*cos(theta2))/sin(theta1);
    
    % 生成直接雷达信道向量hdt (M×1)
    % 包含：路径损耗、阵列相位响应和功率归一化
    hdt = sqrt(10^(-3)*dt^(-alpha_t))*exp(-1j*(0:1:M-1)'*pi*sin(theta1))/sqrt(M);
    
    % 生成RIS-目标信道向量hrt (N×1)
    % 包含：路径损耗和RIS阵列相位响应
    hrt = sqrt(10^(-3)*drt^(-alpha_rt))*exp(-1j*(0:1:N-1)'*pi*sin(theta2));
    
    % Rician衰落模型参数设置
    kappa = 10^(3/10);  % Rician因子，决定LOS(视距)和NLOS(非视距)分量的功率比例
    
    % 生成RIS-基站信道矩阵G的LOS(视距)分量
    % 包含：路径损耗、角度依赖性和阵列响应
    GLos = sqrt(kappa/(1+kappa))*sqrt(10^(-3)*dg^(-alpha_g))*...
           exp(-1j*(0:1:N-1)'*pi*sin(thetar))*...  % RIS阵列相位响应
           exp(-1j*(0:1:M-1)*pi*sin(-thetar))/sqrt(M);  % 基站阵列相位响应与功率归一化
    
    % 生成完整的RIS-基站信道矩阵G (N×M)
    % 结合LOS分量和复高斯分布的NLOS分量
    G = GLos + sqrt(1/(1+kappa))*sqrt(10^(-3)*dg^(-alpha_g))*...
        (randn(N,M)+1i*randn(N,M))/sqrt(2*M);  % 复高斯随机矩阵，功率归一化
    
    % 为每个用户随机生成RIS-用户的入射角，范围在[0, π/2]
    theta_ru = pi/2*rand(K,1);
    
    % 初始化用户相关信道矩阵
    Hu = zeros(K,M);    % 用户到基站的直接信道矩阵 (K×M)
    Hru = zeros(K,N);   % RIS到用户的信道矩阵 (K×N)
    
    % 遍历每个用户，生成相应的信道矩阵
    for k = 1:K
        % 获取当前用户的RIS-用户入射角
        theta_rk = theta_ru(k);
        
        % 计算用户到基站的入射角(基于几何关系)
        theta_k = atan((dg*sin(thetar)-drk*cos(theta_rk))/(dg*cos(thetar)+drk*sin(theta_rk)));
        
        % 计算用户到基站的直接距离
        dk = (dg*sin(thetar)-drk*cos(theta_rk))/sin(theta_k);
        
        % 生成用户k到基站的直接信道向量Hu(k,:) (1×M)
        % 采用Rician衰落模型，包含LOS和NLOS分量
        Hu(k,:) = sqrt(10^(-3)*dk^(-alpha_k))*(...
                  sqrt(kappa/(1+kappa))*exp(-1j*(0:1:M-1)*pi*sin(theta_k))/sqrt(M) + ...  % LOS分量
                  sqrt(1/(1+kappa))*(randn(1,M)+1i*randn(1,M))/sqrt(2*M));  % NLOS分量
        
        % 生成RIS到用户k的信道向量Hru(k,:) (1×N)
        % 同样采用Rician衰落模型，包含LOS和NLOS分量
        Hru(k,:) = sqrt(10^(-3)*drk^(-alpha_rk))*(...
                   sqrt(kappa/(1+kappa))*exp(-1j*(0:1:N-1)*pi*sin(theta_rk)) + ...  % LOS分量
                   sqrt(1/(1+kappa))*(randn(1,N)+1i*randn(1,N))/sqrt(2));  % NLOS分量
    end
    
    % 将所有生成的信道矩阵封装到结构体中返回
    Channel.hdt = hdt;   % 目标-基站直接信道
    Channel.hrt = hrt;   % RIS-目标信道
    Channel.G = G;       % RIS-基站信道
    Channel.Hu = Hu;     % 用户-基站直接信道
    Channel.Hru = Hru;   % RIS-用户信道
end