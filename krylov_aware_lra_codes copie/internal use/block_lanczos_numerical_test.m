clc
clear

%---INPUTS---
% Set up parameters for test
s = 40;

% Set up matrix
n = 4000;
Q = gallery('orthog',n,1);
eigvals = 1:n;
A = Q*diag(eigvals)*Q';
Afun = @(X) A*X;

% Generate block
Omega = randn(n,20);

% Run test
[Q1,T,R,orthogonality_loss,three_term_rr_loss] = block_lanczos(Afun,Omega,s,0,true);

% Plot results
semilogy(orthogonality_loss,'b','LineWidth',3)
hold on
semilogy(three_term_rr_loss,'r','LineWidth',3)
xlabel('Iteration')
ylabel('Frobenius norm error')
legend({'Orthogonality loss','Three term recurrence relation loss'},'Location','best')
hold off
