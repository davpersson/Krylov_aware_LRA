clc
clear

k =10;
eigvals = 0:0.5:(4*pi);
f = @(x) sin(x)+1;
q = 10;

feigvals = f(eigvals);
feigvals_sorted = sort(feigvals,'descend');
flambdakp1 = feigvals_sorted(k+1);
polynomial = chebyshevT(q,feigvals/flambdakp1);

semilogy(eigvals,feigvals,'k')
hold on
semilogy(eigvals,abs(polynomial),'r')
semilogy(eigvals,flambdakp1*ones(length(eigvals),1),'b')
hold off