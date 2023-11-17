function [U,S] = krylov_aware(Afun,f,Omega,s,r,varargin)

l = size(Omega,2);

%Orthogonalize sketch
Omega = orth(Omega);

%Run block Lanczos
[Q1,T] = block_lanczos(Afun,Omega,s,r);

%Compute f(T)
[V,S] = eig(T);
fT = V*diag(f(diag(S)))*V';

%Obtain top (1,1) block of f(T)
b = size(Q1,2);
X = fT(1:b,1:b);
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
