function [ final_est ] = DBMTSpectrogram( y,W,K,U,time_halfbandwidth,num_tapers,Fs )
%[ final_est ] = DBMTSpectrogram( y,W,K,U,time_halfbandwidth,num_tapers,Fs ) computes DBMT
% spectrogram of time-series data. 
%
% Outputs:
%   final_est = DBMT estimate
%
%Inputs:
%   y = Data
%   W = window length in seconds
%   K = # of frequency bin
%   U = # of frequency points to be estimated
%   time_halfbandwidth = N*B
%   num_seq = # of tapers to be used
%   Fs = Sampling Frequency

W = W*Fs;
N = floor(length(y)/W); %# of windows
TOL = 3*10^-3;
final_est = 0;
max_iter = 40;          % usually takes <10 iterations

sigma2 = 0.0001;        % sigma2 and Q^{[0]} needs to be chosed for faster convergence
alpha = 0;

%**********************generate dpss_sequences*****************************
seq_length = W;
[dps_seq,lambda] = dpss(seq_length,time_halfbandwidth,num_tapers);

offset = 1;     
%**************************Fourier matrix**********************************
F = zeros(W,2*U,N);
    for j = 1:N
        thetac = ((j-1)*W+(0:W-1))'*(1:U);
        thetas = ((j-1)*W+(0:W-1))'*(K/2+(K/2-1-U+1:K/2-1));
        % thetac = ((j-1)*W+(0:W-1))'*(offset:U+offset-1);
        % thetas = ((j-1)*W+(0:W-1))'*(K/2+(K/2-U+1-offset:K/2-offset));
        F(:,:,j) = [cos(2*pi/K*thetac), sin(2*pi/K*thetas)]; %/sqrt(U+1)
    end
 
for seq_num = 1:1:num_tapers
    processed_data = Window_then_Taper(y,dps_seq(:,seq_num));
    % x = 0*ones(2*U+1,N);      %k|k or k|N
    x = 0*ones(2*U,N);          %k|k or k|N
    tic
    [x_sol,alpha] = DBMT_EM(seq_num,x,F,N,W,U,processed_data,sigma2,alpha,TOL,max_iter);
    alpha = 0;
    toc
    est = zeros(U,N);
    for k = 1:N
        % est(:,k) = [(x_sol(1,k))^2;((x_sol(2:U+1,k)).^2+flipud((x_sol(U+2:2*U+1,k).^2)))/4];
        est(:,k) = ((x_sol(1:U,k)).^2+flipud((x_sol(U+1:2*U,k).^2)))/4;
    end
    final_est = final_est + lambda(seq_num) * est;
end
final_est = final_est/num_tapers;
%****************************Plot upto 20 Hz*******************************
K = W;
W1 = floor(length(y)/K - W/K);
pcolor(((0:1:W1)*K+W/2)/Fs,Fs*(1:(20*W/Fs))/W,10*log10(W^2*final_est(1:20*W/Fs,:)/Fs));
ylim([0 20]);
shading flat;
colormap('jet');
colorbar;
xlabel('Time(s)','Interpreter','Latex');
ylabel('Frequency(Hz)','Interpreter','Latex');
title('DBMTSpectogram','Interpreter','Latex');
drawnow
end

