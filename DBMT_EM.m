function [ x_sol, alpha_sol ] = DBMT_EM( num,x,F,N,W,U,data,sigma2,alpha,TOL,max_iter )
%[ x_sol, alpha_sol ] = DBMT_EM( num,x,F,N,W,U,data,sigma2,alpha,TOL,max_iter )
% implements EM algorithm required for DBMT analysis.
%
% Outputs:
%   x_sol = time-frequency representation
%   alpha_sol = Estimated Smoothing Parameter
%
% Inputs:
%   num = sequence number
%   data = data
%   x = initial guess
%   F = Fourier matrix
%   N = # of windows
%   W = Window length 
%   U = # of frequency bins
%   sigma2 = \sigma^2, observation noise variance
%   alpha = initial guess for smoothing parameter(usually 0)
%   TOL = tolerence for convergence
%   max_iter = maximum # of of iteration

x0 = 0*ones(2*U,N);    %Dummy 
% x = 0*ones(2*U+1,N);     %k|k or k|N
x1 = 0*ones(2*U,N);    %k|k-1 
K1 = zeros(2*U,W,N);
    
P = zeros(2*U,2*U,N);   %k|k or k|N
for i = 1:N
    % P(:,:,i) = 0.00001*eye(2*U);
    P(:,:,i) = 0.000001*diag(rand(2*U,1)+0.01);
end
P0 = P;             %k|k store
P1 = P;             %k|k-1
P2 = P;             %k-1,k|N
P2(:,:,1)  = zeros(2*U,2*U);
        
Q = zeros(2*U,2*U,N);   %k|k or k|N
for i = 1:N
    Q(:,:,i) = 0.000001*diag(rand(2*U,1)+0.01);
end

B = zeros(2*U,2*U,N);

variable = [];
    
for r = 1:max_iter
%*********************************E step***********************************
    % E step I
    for k = 1:N
        if k == 1
            % one step prediction
            x1(:,1) = alpha*zeros(size(x1(:,1))); 
            % one step variance prediction                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
            P1(:,:,1) =  Q(:,:,1);
        else
            % one step prediction
            x1(:,k) = alpha*x(:,k-1); 
            % one step variance prediction                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
            P1(:,:,k) = alpha^2*P(:,:,k-1) + Q(:,:,k);
        end
        % K1(:,:,k) = P1(:,:,k)*F(:,:,k)'*((F(:,:,k)*P1(:,:,k)*F(:,:,k)'+sigma2*eye(W))\eye(W));
        K1(:,:,k) = ((F(:,:,k)'*F(:,:,k)/sigma2+(P1(:,:,k)\eye(2*U))))\F(:,:,k)'/sigma2;
        % Posteror mode
        x(:,k) = x1(:,k) + K1(:,:,k)*(data(:,k)-F(:,:,k)*x1(:,k));
        % posterior variance
        P(:,:,k) = P1(:,:,k) - K1(:,:,k)*F(:,:,k)*P1(:,:,k);
    end
    P0 = P;
    % E step II FIS
    for k = N-1:-1:1
    % B(:,:,k) = alpha*P(:,:,k)*(P1(:,:,k+1)\eye(2*U+1));
    B(:,:,k) = alpha*P(:,:,k)*(P1(:,:,k+1)\eye(2*U));
    % z(k|N)
    x(:,k) = x(:,k) + B(:,:,k)*(x(:,k+1)-x1(:,k+1));
    % P(k|N)
    P(:,:,k) = P(:,:,k) + B(:,:,k)*(P(:,:,k+1)-P1(:,:,k+1))*(B(:,:,k)'); 
    end
    % E step III state-space covar algorithm
    for k = 2:N
        P2(:,:,k) = B(:,:,k-1)*P(:,:,k);
    end
%*********************************M step***********************************
%****************************Update alpha**********************************
    W1 = 0;
    for k = 1:N-1
        % W1 = W1 + x(:,k)'*(Q(:,:,k+1)\x(:,k)) + trace(Q(:,:,k+1)\P(:,:,k));
        W1 = W1 + x(:,k).^2./diag(Q(:,:,k+1)) + diag(P(:,:,k))./diag(Q(:,:,k+1));
    end
    W1 = sum(W1);
    W2 = 0;
    for k = 2:N
        % W2 = W2 + x(:,k-1)'*(Q(:,:,k)\x(:,k)) + trace(Q(:,:,k)\(P2(:,:,k)+P2(:,:,k)')/2);
        W2 = W2 + (x(:,k-1).*x(:,k))./diag(Q(:,:,k)) + diag(P2(:,:,k))./diag(Q(:,:,k)); 
    end
    W2 = sum(W2);
    alpha = W2/W1;


%*************************Stopping condition*******************************
    variable = [variable, norm(x-x0,'fro')/norm(x,'fro')];
    if variable(r) < TOL
        break;      % Convergence reached
    else
%***********************GOTO next iterations*******************************
        x0 = x;
        P = P0;
        Q(:,:,1) = diag(diag(x(:,1)*x(:,1)'+P(:,:,1)));
        P2(:,:,1) = Q(:,:,1);
        for k = 2:N
            % W3 = P2(:,:,k) + x(:,k)*x(:,k-1)';
            Q(:,:,k) = diag((x(:,k)-alpha*x(:,k-1)).^2 + diag(P(:,:,k))+alpha^2*diag(P(:,:,k-1))-2*alpha*diag(P2(:,:,k)));
            % Q(:,:,k) = ((x(:,k)*x(:,k)'+P(:,:,k) + alpha^2*(x(:,k-1)*x(:,k-1)'+P(:,:,k-1))-alpha*(W3+W3')));
        end
    end
end

x_sol = x;
alpha_sol = alpha;

%*********************Save Variables for CI********************************
pathname = fileparts('./DBMT/file');
file_name = sprintf('taper%d.mat',num);
matfile = fullfile(pathname, file_name);
save(matfile, 'x', 'P');

end

