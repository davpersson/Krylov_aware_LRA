function Y = quadform(Afun,f,Omega,q)

b = size(Omega,2);
[Q,T,R1] = block_lanczos(Afun,Omega,q,0);
[U,S] = eig(T);
fT = U*diag(f(diag(S)))*U';
Y = R1'*fT(1:b,1:b)*R1;


end