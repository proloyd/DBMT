function [ x_sol, alpha_sol ] = log_DBMT_EM( num,x,N,Y,alpha,TOL,max_iter )
%EM_algo Computes time-frequency representation given the data, Y
%   num = sequence number
%   Y = data
%   x = initial guess
%   N = # of windows
%   alpha = smoothing parameter
%   TOL = tolerence for convergence
%   max_iter = maximum # of of iteration

x0 = 0*ones(size(Y));    %Dummy 
% x = 0*ones(2*U+1,N);     %k|k or k|N
x1 = 0*ones(size(Y));    %k|k-1 
    
P = zeros(size(Y,1),size(Y,1),N);   %k|k or k|N
for i = 1:N
    P(:,:,i) = 10^-1*eye(size(Y,1));
end
P0 = P;             %k|k store
P1 = P;             %k|k-1
P2 = P;             %k-1,k|N
P2(:,:,1)  = zeros(size(Y,1),size(Y,1));
    
Q = zeros(size(Y,1),size(Y,1),N);   %k|k or k|N
for i = 1:N
    Q(:,:,i) = 10^-1*eye(size(Y,1));
end

B = zeros(size(Y,1),size(Y,1),N);

variable = [];

nu = 2; % initialize nu
    
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
        % Posteror mode and posterior variance
        [x(:,k),P(:,:,k)] = Post_mode_var(Y(:,k), x1(:,k), nu, P1(:,:,k));
    end
    P0 = P;
    % E step II FIS
    for k = N-1:-1:1
    B(:,:,k) = alpha*P(:,:,k)*(P1(:,:,k+1)\eye(size(Y,1)));
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
    
%*****************************Update nu************************************
    cnst = 0;
    for k = 1:N
        cnst = cnst + sum(Y(:,k)-x(:,k)-exp(Y(:,k)-x(:,k)));
    end
    cnst = cnst/(N*size(Y,1));
    f = @(nu) log(gamma(nu)) - nu*log(2*nu) - (cnst+1)*nu;
    df =@(nu) psi(nu) - log(nu) - (cnst+1);
    d2f =@(nu) psi(1,nu) - 1/nu;
    tau = 1;
    nu = 1;
    dir = -df(nu)/d2f(nu);
    for k = 1:100
        while  nu+tau*dir < 0.01 
            tau = tau/2;
            if tau < 10^-10
                break
            end
        end
        nu = nu + tau * dir;
        dir = -df(nu)/d2f(nu);
        if -df(nu)*dir < 10^-4  || tau < 10^-10
            break
        end
        tau = 1;
    end
    fprintf('itr = %d, nu = %f \n', r,nu);
          

%*************************Stopping condition*******************************
    variable = [variable, norm(x-x0,'fro')/norm(x,'fro')];
    if variable(end) < TOL
        break;      % Convergence reached
    else
%***********************GOTO next iterations*******************************
        x0 = x;
%***********************************Update Q*******************************
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
pathname = fileparts('./log_DBMT/file');
file_name = sprintf('taper%d.mat',num);
matfile = fullfile(pathname, file_name);
save(matfile, 'x', 'P');

end

