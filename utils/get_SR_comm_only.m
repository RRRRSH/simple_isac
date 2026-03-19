function [W,phi,Vsr,gammat_r] = get_SR_comm_only(Prms,Channel,phi,W)
% get_SR_comm_only - 无雷达约束的纯通信优化（获取系统的 Sum-Rate 物理上限）
% 本函数是为 Pareto 边界测试提供的终极天花板参考值

M = Prms.M; N = Prms.N; K = Prms.K; sigmar2 = Prms.sigmar2; sigmak2 = Prms.sigmak2;
sigmat2 = Prms.sigmat2; Nmax = Prms.Nmax; res_th = Prms.res_th; nL = Prms.L;
P = Prms.P; hdt = Channel.hdt; hrt = Channel.hrt; G = Channel.G; Hu = Channel.Hu; Hru = Channel.Hru;

Hk = Hu + Hru*diag(phi)*G;
r = zeros(K,1); c = zeros(K,1);
for k = 1:K
    r(k) = abs(Hk(k,:)*W(:,k))^2/(norm(Hk(k,:)*W,2)^2-abs(Hk(k,:)*W(:,k))^2+sigmak2);
    c(k) = sqrt(1+r(k))*Hk(k,:)*W(:,k)/(norm(Hk(k,:)*W,2)^2+sigmak2);
end

iter = 1; res = 1;
Vsr = zeros(1,Nmax);

while iter <= Nmax && res >= res_th
    fprintf('        -> [Comm-only] 迭代 %d: 优化 RIS 相位... ', iter);
    D = zeros(N,N); g = zeros(N,1);
    for k = 1:1:K
        g = g + 2*sqrt(1+r(k))*c(k)*diag(W(:,k)'*G')*Hru(k,:)';
        for j = 1:1:K
            temp = diag(W(:,j)'*G')*Hru(k,:)';
            D = D + abs(c(k))^2*(temp*temp');
            g = g - 2*abs(c(k))^2*diag(W(:,j)'*G')*Hru(k,:)'*Hu(k,:)*W(:,j);
        end
    end
    
    phi = get_initial_phi(-D/norm(g),g/norm(g));
    fprintf('完成. 优化 基站波束... ');

    Hk = Hu + Hru*diag(phi)*G;
    a = zeros(M*(K),1); B = zeros(K*(K),M*(K));
    for k = 1:1:K
        a((k-1)*M+1:k*M) = 2*sqrt(1+r(k))*c(k)*Hk(k,:)';
        for j = 1:1:K
            Tj = zeros(M,M*(K)); Tj(:,(j-1)*M+1:j*M) = eye(M);
            B((k-1)*(K)+j,:) = abs(c(k))*Tj.'*Hk(k,:).';
        end
    end
    
    cvx_begin quiet
    cvx_solver SeDuMi
    variable w(M*K,1) complex
    minimize square_pos(norm(B*w,2))-real(a'*w)
    subject to
    norm(w,2) <= sqrt(P);
    cvx_end
    
    if ~strcmpi(cvx_status, 'Solved') && ~strcmpi(cvx_status, 'Inaccurate/Solved')
        fprintf('无解! 强制退出.\n'); break;
    end
    fprintf('完成.\n');

    W = reshape(w,M,K);
    r = zeros(K,1); c = zeros(K,1);
    for k = 1:1:K
        r(k) = abs(Hk(k,:)*W(:,k))^2/(norm(Hk(k,:)*W,2)^2-abs(Hk(k,:)*W(:,k))^2+sigmak2);
        c(k) = sqrt(1+r(k))*Hk(k,:)*W(:,k)/(norm(Hk(k,:)*W,2)^2+sigmak2);
    end

    Vsr(iter) = sum(log2(1+r));
    if iter > 10
        res = (Vsr(iter)-Vsr(iter-1))/Vsr(iter-1);
    end
    iter = iter + 1;
end

Vsr(iter:end) = [];
Ht = (hdt + G.'*diag(phi)*hrt)*(hdt.' + hrt.'*diag(phi)*G);
u = kron(eye(K),Ht)*vec(W)/((vec(W))'*(kron(eye(K),Ht'*Ht))*vec(W));
gammat_r = real(nL*sigmat2*abs(u'*kron(eye(K),Ht)*vec(W))^2/(sigmar2*u'*u));

fprintf('      > Comm-only 已收敛. 终态: SR=%.4f, 自然溢出SNR=%.2f dB\n', Vsr(end), 10*log10(abs(gammat_r)+1e-12));
end
