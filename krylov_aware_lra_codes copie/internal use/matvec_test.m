clc
clear

%---INPUTS---
% Set up parameters for test
s_list = 5:20;

% Set up matrix function
%beta = 0.5;
%f = @(x) exp(-beta*x);
f = @(x) exp(x);
%f = @(x) x.^4;

% Set up matrix
% n = 4000;
% Q = gallery('orthog',n,1);
% eigvals = 1:n;
% eigvals = abs(randn(n,1));
% A = Q*diag(eigvals)*Q';
% fA = Q*diag(f(eigvals))*Q';
% Afun = @(X) A*X;
% feigvals = sort(f(eigvals),'descend');

A = sparse(uq_laplaceeig(0.01,1));
n = size(A,1);
[U,S] = eig(full(A));
fA = U*diag(f(diag(S)))*U';
Afun = @(X) A*X;
fAfun_exact = @(X) fA*X;

Omega = randn(n,70);
Yexact = fA*Omega;

error = zeros(1,length(s_list));


iteration = 0;
for s = s_list
    
    iteration = iteration + 1
    Y = matvec(Afun,f,Omega,s);
    error(iteration) = norm(Yexact - Y,'fro')/norm(Yexact,'fro');
    
    
end

semilogy(s_list,error,'k-*','LineWidth',3)
xlabel('Iterations')
ylabel('Relative Frobenius norm error')
set(gca,'Fontsize',14)
