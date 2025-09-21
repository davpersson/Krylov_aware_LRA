% Run experiments to produce results for Section 5.3

clc
clear
rng(0)

%-- Add paths --
addpath('results')
addpath('codes')

%-- Set parameters --
% Determine what to keep constant in the two methods
constant = 'basis'; % The single vector method produces a basis that has the same dimension as the block method. Extra effort spent on approximating quadratic form. 
%constant = 'quadratic_form'; % The single vector method produces a basis that has a larger dimension than the block method. No extra effort spent on approximating quadratic form.

tic
for matrix = 1:4
    
    matrix
    
    %-- Select matrix --
    if matrix == 1
       
        k = 60;
        s_list = 38:42;
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
        s_list = 2:20;
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
        k = 10
        s_list = 2:30;
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
    svk_krylov_aware_untruncated = zeros(1,length(s_list));
    svk_krylov_aware_truncated = zeros(1,length(s_list));
    
    % Check that there is a 1-to-1 correspondence between r_list and s_list
    if length(s_list)~=length(r_list)
    
        error('The lengths of s_list and r_list must be the same')
    
    end
    
    % Generate random matrix
    Omega = randn(n,k);
    omega_single = randn(n,1);
    ell = size(Omega,2);
    
    
    for i = 1:length(s_list)

        fprintf('%i / %i \n',i,length(s_list))
        
        % Number of iterations
        s = s_list(i);
        r = r_list(i);
        
        % Krylov aware error
        [U,S] = krylov_aware(Afun,f,Omega,s,r);
        krylov_aware_untruncated(i) = norm(fA - U*S*U','fro')/normfA;
        krylov_aware_truncated(i) = norm(fA - U(:,1:k)*S(1:k,1:k)*U(:,1:k)','fro')/normfA;
        
        % single vector Krylov aware error
        if strcmp(constant,'basis')
            s_hat = ell*s;
            r_hat = ell*r;
        elseif strcmp(constant,'quadratic_form')
            s_hat = ell*s+(ell-1)*r;
            r_hat = r;
        end
        [U,S] = svk_krylov_aware(Afun,f,omega_single,s_hat,r_hat);
        svk_krylov_aware_untruncated(i) = norm(fA - U*S*U','fro')/normfA;
        svk_krylov_aware_truncated(i) = norm(fA - U(:,1:k)*S(1:k,1:k)*U(:,1:k)','fro')/normfA;
        
        
    end

    filename = append(filename,'_svk',constant);

    save(filename,'s_list','r_list','optimal','krylov_aware_untruncated',...
        'krylov_aware_truncated','svk_krylov_aware_truncated','svk_krylov_aware_untruncated','k','ell')
    
    %-- Plot results --
    svk_error_plotter(filename)
    
end
toc