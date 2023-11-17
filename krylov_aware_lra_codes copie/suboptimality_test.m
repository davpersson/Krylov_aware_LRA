clc
clear
rng(0)

%-- Add paths --
addpath('results')
addpath('codes')

%-- Set parameters --
p = 0;
qmax = 10;

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
        optfA = norm(feigvals(k+1:end));
        filename = 'results/exponential_integrator_suboptimality';
        
        
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
        optfA = norm(feigvals(k+1:end));
        filename = 'results/estrada_suboptimality';
        
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
        Afun = @(X) A*X;
        fAfun_exact = @(X) fA*X;
        optfA = norm(feigvals(k+1:end));
        filename = 'results/quantum_spin_suboptimality';
        
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
        Afun = @(X) A*X;
        fAfun_exact = @(X) fA*X;
        optfA = norm(feigvals(k+1:end));
        filename = 'results/synthetic_log_suoboptimality';
        
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
        Afun = @(X) A*X;
        fAfun_exact = @(X) fA*X;
        optfA = norm(feigvals(k+1:end));
        filename = 'other_suboptimality';
        
    end
    
    %-- Run test --
    krylov_aware_ = zeros(1,length(s_list));
    subspace_iteration = zeros(1,length(s_list));
    q_optimal = zeros(1,length(s_list));
    
    % Check that there is a 1-to-1 correspondence between r_list and s_list
    if length(s_list)~=length(r_list)
    
        error('The lengths of s_list and r_list must be the same')
    
    end
    
    % Generate random matrix
    Omega = randn(n,k+p);
    
    for i = 1:length(s_list)
        
        % Number of iterations
        s = s_list(i);
        r = r_list(i);
        
        % Krylov aware untruncated error
        [U,S] = krylov_aware(Afun,f,Omega,s,r);
        krylov_aware_(i) = norm(fA - U(:,1:k)*S(1:k,1:k)*U(:,1:k)','fro')/optfA - 1;
        
        epsilon_opt = inf;
        qopt = 1;
        
        b = size(Omega,2);

        %Run block Lanczos
        [Q1,T,R1] = block_lanczos(Afun,Omega,s,0);
        [V,D] = eig(T);
        
        for q = 1:qmax
           

            %Compute f(T)
            fT = V*diag(f(diag(D)).^q)*V';
            Y = Q1*fT(:,1:b)*R1;
            Q = orth(Y);
            
            [U,S] = eig(Q'*matvec(Afun,f,Q,r));
            U = Q*U;
            
            % Sort eigenvalues
            [~,ind] = sort(abs(diag(S)),'descend');
            U = U(:,ind);
            S = S(ind,ind);
            
            k_ = size(U,2);
            k__ = min(k,k_);
            epsilon = norm(fA - U(:,1:k__)*S(1:k__,1:k__)*U(:,1:k__)','fro')/optfA - 1;
            
            if epsilon <= epsilon_opt
                
                epsilon_opt = epsilon;
                qopt = q;
                
            end
            
            subspace_iteration(i) = epsilon_opt;
            q_optimal(i) = qopt;
        end
        
    end
    
    save(filename,'s_list','r_list','krylov_aware_',...
        'subspace_iteration','q_optimal','k')
    
    %-- Plot results --
    suboptimality_plotter(filename)
    
end
toc