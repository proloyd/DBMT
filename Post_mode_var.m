function [ x_sol, Sigma_sol ] = Post_mode_var( y, x, nu, Sigma )
% Finds the posterior mode and variance
%   y = observation
%   x = x_{(k|k-1)} One step prediction
%   Sigma = \Sigma_{k|k-1} Prediction variance

MAXLEN = 100;
TOL = 10^-6;

% degrees of freedom
% nu = 1.0;

% function def
f = @(x_sol) nu*sum(x_sol-y) + sum(exp(y-x_sol+log(nu))) + (x_sol-x)'*(Sigma\(x_sol -x));
grad = @(x_sol) nu*ones(size(y)) - exp(y-x_sol+log(nu)) + Sigma\(x_sol-x);
Hes = @(x_sol) Sigma\eye(length(y)) + diag(exp(y-x_sol+log(nu)));

% f = @(x_sol) nu*sum(x_sol-y) + sum(exp(y+log(nu)).*(exp(-x_sol)+0.5*exp(-x_sol).*diag(Sigma))) + (x_sol-x)'*(Sigma\(x_sol -x));
% grad = @(x_sol) nu*ones(size(y)) - exp(y+log(nu)).*(exp(-x_sol)+0.5*exp(-x_sol).*diag(Sigma)) + Sigma\(x_sol-x);
% Hes = @(x_sol) Sigma\eye(length(y)) + diag(exp(y+log(nu)).*(exp(-x_sol)+0.5*exp(-x_sol).*diag(Sigma)));


x_sol = zeros(size(y));
g = grad(x_sol);
H = Hes(x_sol);
d = H\g;
tau = 1;  % Stepsize
alpha = 0.01;

res = g'*d;

for i = 1:MAXLEN
    % Backtracking Line search
    while f(x_sol - tau*d) > f(x_sol) - alpha * tau * (g'*d)
        tau = tau/2;
    end
    x_sol = x_sol - tau*d;
    g = grad(x_sol);
    H = Hes(x_sol);
    d = H\g;
    if tau < 10^-15
        break;
    end
    tau = 1;  % Stepsize
    res = [res; g'*d];
    if res(end)/res(1) < TOL
        break;
    end
end
Sigma_sol = Hes(x_sol)\eye(length(y));

end

