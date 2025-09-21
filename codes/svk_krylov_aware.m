function [U,S] = svk_krylov_aware(Afun,f,omega,s,r,varargin)
% Code for single vector Krylov aware low rank approximation

%Orthogonalize sketch
omega = omega/norm(omega);

%Run Lanczos
%[Q1,T] = lanczos(Afun,omega,s,r);
A = Afun(sparse(eye(size(omega,1))));
[Q1,T] = lanczos_through_rktoolbox(A,omega,s,r);

T = (T + T')/2;

%Compute f(T)
[V,S] = eig(T);
fT = V*diag(f(diag(S)))*V';

%Obtain top (1,1) block of f(T)
b = size(Q1,2);
X = fT(1:b,1:b); X = (X+X')/2;
[V,S] = eig(X);
U = Q1*V;

% Sort eigenvalues
[~,ind] = sort(abs(diag(S)),'descend');
U = U(:,ind);
S = S(ind,ind);

%Truncate
if nargin == 6
    k = varargin{1};
    U = U(:,1:k);
    S = S(1:k,1:k);
end

end
