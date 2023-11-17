function [Q_stack,T,R1] = block_lanczos(Afun,Z,q)

b = size(Z,2);

[Q_current,R1] = qr(Z,0);
R = R1;
Q_stack = [];

for k = 1:q
    
    if k == 1
        
        Z = Afun(Q_current);
        
    else
        
        Z = Afun(Q_current) - Q_prev*R;
        
    end
    
    M = Q_current'*Z;
    
    if k == 1
        
        T = M;
        
    else
       
        m = size(T,1);
        T = [T zeros(m,b);zeros(b,m) M];
        b_old = size(R,2);
        T((end-b+1):end,(end-b-b_old+1):(end-b)) = R;
        T((end-b-b_old+1):(end-b),(end-b+1):end) = R';
        
    end
    
    if k == q
        
        Q_stack = [Q_stack Q_current];
        return
        
    end
    
    Z = Z - Q_current*M;

    if k > 1
        Z = Z - Q_stack*(Q_stack'*Z);
    end
    
    Q_prev = Q_current;
    [Q_current,R,P] = qr(Z,0);
    invP = zeros(1,length(P));
    invP(P) = 1:length(P);
    R = R(:,invP);
    r = rank(R);
    Q_stack = [Q_stack Q_prev];
    
    if norm(R,'fro') < 1e-10

        return
        
    end

    if r < b
        
        R = R(1:r,:);
        Q_current = Q_current(:,1:r);
        b = r;
        

    end
    
end

end