function [Q1,T,R0,varargout] = block_lanczos(Afun,Z,s,r,varargin)

if nargin == 5 && varargin{1} == true
    
    orthogonality_loss = zeros(1,s);
    three_term_rr_loss = zeros(1,s);
    
end

q = s + r;

% Find block size
b = size(Z,2);

% First iteration
[Q,R0] = qr(Z,0);
R = R0;

for k = 0:q-1
    
    %Set Qk
    Qk = Q(:,(1+k*b):((k+1)*b));
    
    if k == 0
        
        Z = Afun(Qk);
        
    else
        
        % Set Qk-1
        Qkm1 = Q(:,(1+(k-1)*b):(k*b));
        
        Z = Afun(Qk) - Qkm1*R';
        
    end
    
    % Obtain diagonal block in T
    M = Qk'*Z;
    
    if k == 0
        
        T = M;
        
    else
        
        m = size(T,1);
        T = [T zeros(m,b);zeros(b,m) M];
        b_old = size(R,2);
        T((end-b+1):end,(end-b-b_old+1):(end-b)) = R;
        T((end-b-b_old+1):(end-b),(end-b+1):end) = R';
        
    end
    
    if k <= s-1
        
        Q1 = Q;
        
    end
    
    if k == q-1
        
        if nargin == 5 && varargin{1} == true
            varargout{1} = orthogonality_loss;
            varargout{2} = three_term_rr_loss;
        end
        return
        
    end
        
    
    % Reorthogonalization
    Z = Z - Qk*M;
    
    % Double reorthgonalization
    if k > 0
        
        Z = Z - Q(:,1:(k*b))*(Q(:,1:(k*b))'*Z);
        
    end
    
    % Obtain next block
    [Qkp1,R,P] = qr(Z,0);
    
    % Check if rank deficient
    if min(abs(diag(R))) <= (1e-12)*max(abs(diag(R)))
        
        % New block size
        b = max(find(abs(diag(R)) > (1e-10)*max(abs(diag(R)))));
        
        % Truncate
        R = R(:,1:b);
        Qkp1 = Qkp1(:,1:b);
        
    end
    
    % Permute back
    invP = zeros(1,length(P));
    invP(P) = 1:length(P);
    R = R(:,invP);
    if nargin == 5 && varargin{1} == true && k <= s-1
        
        orthogonality_loss(k+1) = norm(Q'*Q-eye(size(Q,2)),'fro');
        Ek = eye(size(Q,2)); 
        Ek = Ek(:,((end-b)+1):end);
        three_term_rr_loss(k+1) = norm(Afun(Q) - Q*T-Qkp1*R*Ek');
        
    end
        
        
    Q = [Q Qkp1];
    
    
end

end