function [U,D] = randSVD(Afun,Omega,varargin)

% Randomized SVD
Q = orth(Afun(Omega));

[U,D] = eig(Q'*Afun(Q));
U = Q*U;

% Sort eigenvalues
[~,ind] = sort(abs(diag(D)),'descend');
U = U(:,ind);
D = D(ind,ind);

% Truncate
if nargin == 4
    k = varargin{1};
    U = U(:,1:k);
    D = D(1:k,1:k)
end
end