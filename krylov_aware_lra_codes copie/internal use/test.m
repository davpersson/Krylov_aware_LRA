clc
clear

%---INPUTS---
% Set up parameters for test
k_list = [10,20,30,40];
p = 4;
s_list = 5:20;
r_list = s_list;

% Set up matrix function
beta = 1;
f = @(x) exp(-beta*x);

% Set up matrix
n = 4000;
Q = gallery('orthog',n,1);
eigvals = 2*log(1:n);
A = Q*diag(eigvals)*Q';
fA = Q*diag(f(eigvals))*Q';
Afun = @(X) A*X;
feigvals = sort(f(eigvals),'descend');

%---TESTS---
% Allocate space to save results
optimal = zeros(length(k_list),length(s_list));
krylov_aware_untruncated = zeros(length(k_list),length(s_list));
krylov_aware_truncated = zeros(length(k_list),length(s_list));
randSVD_untruncated = zeros(length(k_list),length(s_list));
randSVD_truncated = zeros(length(k_list),length(s_list));

% Check that there is a 1-to-1 correspondence between r_list and s_list
if length(s_list)~=length(r_list)
    
    error('The lengths of s_list and r_list must be the same')
    
end

for i = 1:length(k_list)
    
    % Rank
    k = k_list(i)
    
    % Compute optimal low rank approximation error
    optimal(i,:) = norm(feigvals((k+1):end))*ones(1,length(s_list));
    
    % Generate random matrix
    Omega = randn(n,k+p);
    
    for j = 1:length(s_list)
        
        % Number of iterations
        s = s_list(j);
        r = r_list(j);
        
        % Krylov aware low rank approximation error
        [U,S] = krylov_aware(Afun,f,Omega,s,r);
        krylov_aware_untruncated(i,j) = norm(fA - U*S*U','fro');
        krylov_aware_truncated(i,j) = norm(fA - U(:,1:k)*S(1:k,1:k)*U(:,1:k)','fro');
        
        % Krylov aware low rank approximation error
        [U,S] = krylov_aware(Afun,f,Omega,s+r,0);
        krylov_aware_untruncated2(i,j) = norm(fA - U*S*U','fro');
        krylov_aware_truncated2(i,j) = norm(fA - U(:,1:k)*S(1:k,1:k)*U(:,1:k)','fro');
        
        % randSVD low rank approximation error
        fAfun = @(X) matvec(Afun,f,X,s);
        [U,S] = randSVD(fAfun,Omega);
        randSVD_untruncated(i,j) = norm(fA - U*S*U','fro');
        k_ = size(U,2);
        k__ = min(k,k_);
        randSVD_truncated(i,j) = norm(fA - U(:,1:k__)*S(1:k__,1:k__)*U(:,1:k__)','fro');
        
    end
    
end

%---PLOTS---
for i = 1:length(k_list)
    
    figure(i)
    normfA = norm(fA,'fro');
    semilogy(s_list + r_list,optimal(i,:)/normfA,'k','LineWidth',3)
    hold on
    semilogy(s_list + r_list,krylov_aware_untruncated(i,:)/normfA,'b--*','LineWidth',3)
    semilogy(s_list + r_list,krylov_aware_truncated(i,:)/normfA,'b-*','LineWidth',3)
    semilogy(s_list + r_list,krylov_aware_untruncated2(i,:)/normfA,'g--*','LineWidth',3)
    semilogy(s_list + r_list,krylov_aware_truncated2(i,:)/normfA,'g-*','LineWidth',3)
    semilogy(2*s_list,randSVD_untruncated(i,:)/normfA,'r--*','LineWidth',3)
    semilogy(2*s_list,randSVD_truncated(i,:)/normfA,'r-*','LineWidth',3)
    legend({'Optimal','Krylov aware (untruncated)','Krylov aware (truncated)',...
        'Krylov aware (untruncated, r=0)','Krylov aware (truncated, r=0)',...
        'randSVD (untruncated)','randSVD (truncated)'},'location','best')
    xlabel('Number of matrix vector products with $A$','Interpreter','latex')
    ylabel('Relative Frobenius norm error','Interpreter','latex')
    title_text = append('$k = $',num2str(k_list(i)));
    title(title_text,'interpreter','latex')
    set(gca,'Fontsize',14)
    hold off
    
end
