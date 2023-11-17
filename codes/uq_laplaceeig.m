function L = uq_laplaceeig(kappa,lambda)

%Function to create A_pde

N = 100;

h = 1/N;

T = diag(-2*ones(N-1,1)) + diag(ones(N-2,1),-1) + diag(ones(N-2,1),1);
J = zeros(N,N); J(end,end) = 1/2;
That = T - 2*eye(size(T));
Ttilde = diag(-2*ones(N,1)) + diag(ones(N-1,1),-1) + diag(ones(N-1,1),1);
L = kappa*(kron(eye(N),T) + kron(Ttilde,eye(N-1)) - kron(J,That))/h^2 + ...
    lambda*(eye(N*(N-1)) - kron(J,eye(N-1)));

end