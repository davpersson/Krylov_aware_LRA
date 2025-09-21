function Y = matvec(Afun,f,Omega,q)
% Code for approximating matvecs with f(A) using the block Lanczos method

b = size(Omega,2);
[Q,T,R1] = block_lanczos(Afun,Omega,q,0);
[U,S] = eig(T);
fT = U*diag(f(diag(S)))*U';
Y = Q*fT(:,1:b)*R1;


end