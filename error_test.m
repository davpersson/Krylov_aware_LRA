clc
clear
rng(0)

%-- Add paths --
addpath('results')
addpath('codes')

%-- Set parameters --
% Set the oversampling parameter
%p = 0;
p = 5;

tic
for matrix = 1:4
    
    matrix
    
    %-- Select matrix --
    if matrix == 1
       
        k = 60;
        s_list = 30:50;
        r_list = s_list;
        f = @(x) exp(x);
        A = sparse(uq_laplaceeig(0.01,1));
        n = size(A,1);
        [U,S] = eig(full(A));
        fA = U*diag(f(diag(S)))*U';
        Afun = @(X) A*X;
        fAfun_exact = @(X) fA*X;
        eigvals = sort(diag(S),'descend');
        feigvals = f(eigvals);
        normfA = norm(feigvals);
        optimal = sort(abs(feigvals),'descend'); 
        optimal = norm(optimal(k+1:end))/normfA;
        filename = 'results/exponential_integrator';
        
        
    elseif matrix == 2
        
        k = 10;
        s_list = 5:20;
        r_list = s_list;
        f = @(x) exp(x);
        A = sparse(create_roget_mat());
        n = size(A,1);
        [U,S] = eig(full(A));
        fA = U*diag(f(diag(S)))*U';
        Afun = @(X) A*X;
        fAfun_exact = @(X) fA*X;
        eigvals = sort(diag(S),'descend');
        feigvals = f(eigvals);
        normfA = norm(feigvals);
        optimal = sort(abs(feigvals),'descend'); 
        optimal = norm(optimal(k+1:end))/normfA;
        filename = 'results/estrada';
        
    elseif matrix == 3
        
        N = 14;
        h = 10;
        beta = 0.3;
        
        n = 2^N;
        k = 10;
        s_list = 10:20;
        r_list = s_list;
        f = @(x) exp(-beta*x);
        eigvals = sort(tfim_eigs(N,h) + (1+h)*N,'ascend');
        feigvals = f(eigvals);
        A = diag(sparse(eigvals));
        fA = diag(sparse(feigvals));
        normfA = norm(feigvals);
        Afun = @(X) A*X;
        fAfun_exact = @(X) fA*X;
        optimal = norm(feigvals(k+1:end))/normfA;
        filename = 'results/quantum_spin';
        
    elseif matrix == 4
        
        n = 5000;
        k = 30;
        s_list = 1:5;
        r_list = s_list;
        f = @(x) log(x);
        eigvals = exp((1:n).^(-2));
        feigvals = (1:n).^(-2);
        A = sparse(diag(eigvals));
        fA = sparse(diag(feigvals));
        normfA = norm(feigvals);
        Afun = @(X) A*X;
        fAfun_exact = @(X) fA*X;
        optimal = norm(feigvals(k+1:end))/normfA;
        filename = 'results/synthetic_log';
        
    else
        
        n = 5000;
        k = 30;
        s_list = 1:5;
        r_list = s_list;
        f = @(x) exp(-x);
        eigvals = exp((1:n).^(-2));
        feigvals = (1:n).^(-2);
        A = sparse(diag(eigvals));
        fA = sparse(diag(feigvals));
        normfA = norm(feigvals);
        Afun = @(X) A*X;
        fAfun_exact = @(X) fA*X;
        optimal = norm(feigvals(k+1:end))/normfA;
        filename = 'other';
        
    end
    
    %-- Run test --
    krylov_aware_untruncated = zeros(1,length(s_list));
    krylov_aware_truncated = zeros(1,length(s_list));
    randSVD_truncated = zeros(1,length(s_list));
    %exactrandSVD_truncated = zeros(1,length(s_list));
    %power_iteration_truncated = zeros(1,length(s_list));
    %q_optimal = zeros(1,length(s_list));
    
    % Check that there is a 1-to-1 correspondence between r_list and s_list
    if length(s_list)~=length(r_list)
    
        error('The lengths of s_list and r_list must be the same')
    
    end
    
    % Generate random matrix
    Omega = randn(n,k+p);
    [U,S] = randSVD(fAfun_exact,Omega);
    k_ = size(U,2);
    k__ = min(k,k_);
    exactrandSVD_truncated = norm(fA - U(:,1:k__)*S(1:k__,1:k__)*U(:,1:k__)','fro')/norm(fA,'fro');
    
    
    for i = 1:length(s_list)
        
        % Number of iterations
        s = s_list(i);
        r = r_list(i);
        
        % Set a function to compute f(A)*x
        fAfun = @(X) matvec(Afun,f,X,s);
        
        % Krylov aware untruncated error
        [U,S] = krylov_aware(Afun,f,Omega,s,r);
        krylov_aware_untruncated(i) = norm(fA - U*S*U','fro')/normfA;
        krylov_aware_truncated(i) = norm(fA - U(:,1:k)*S(1:k,1:k)*U(:,1:k)','fro')/normfA;
        
        % randSVD low rank approximation error
        [U,S] = randSVD(fAfun,Omega);
        k_ = size(U,2);
        k__ = min(k,k_);
        %randSVD_untruncated(j) = norm(fA - U*S*U','fro');
        randSVD_truncated(i) = norm(fA - U(:,1:k__)*S(1:k__,1:k__)*U(:,1:k__)','fro')/normfA;
        
        % % Subspace iteration with optimal q
        %opterr = randSVD_truncated(j);
        %qopt = 1;
        %for q = 2:qmax
            
        %    err = 
            
        %end
        
    end
    
    save(filename,'s_list','r_list','optimal','krylov_aware_untruncated',...
        'krylov_aware_truncated','randSVD_truncated','exactrandSVD_truncated','k')
    
    %-- Plot results --
    error_plotter(filename)
    
end
toc