function [y] = Window_then_Taper( x,seq )
%[y] = Window_then_Taper( x,seq ) windows the data
% and then multiply by taper (Used in DBMT algo). 
%
% Output:
%  y = data
%
% Input:
%  x = time-series data
%  seq = Taper


N = length(seq);
W = floor(length(x)/N);

y = zeros(N,W);
for k = 0:1:W-1
   d = x(1+k*N:(k+1)*N); 
   y(:,k+1) = d.*seq;
end

end

