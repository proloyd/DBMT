function [ final_est ] = log_DBMTSpectrogram( y,W,U,time_halfbandwidth,num_tapers,Fs )
%Computed spectogram using the proposed method
%   y = Times series data
%   W = window length in seconds
%   U = # of frequency points to be estimated
%   K = # of frequency bin
%   Fs = sampling frequency

N = floor(length(y)/(W*Fs)); %# of windows
TOL = 3*10^-2;
max_iter = 40;

final_est = 0;

for seq_num = 1:1:num_tapers
    tic
    S = TSpectrogram(y,W,0,time_halfbandwidth,seq_num,Fs);
    Y = log(S(1:U+1,:));     
    x = 0*ones(size(Y));     %k|k or k|N
    alpha = 0;
    est = zeros(size(Y));
    [est,alpha] = log_DBMT_EM(seq_num,x,N,Y,alpha,TOL,max_iter); 
    toc
    final_est = final_est + exp(est);
end
final_est = final_est/num_tapers;

%****************************Plot upto 20 Hz*******************************
W = floor(Fs*W);
K = W;
W1 = floor(length(y)/K - W/K);
pcolor(((0:1:W1)*K+W/2)/Fs,Fs*(0:(20*W/Fs)-1)/N,10*log10(final_est(1:20*W/Fs,:)/Fs));
shading flat;
colormap('jet');
colorbar;
xlabel('Time(s)','Interpreter','Latex');
ylabel('Frequency(Hz)','Interpreter','Latex');
title('log-DBMTSpectogram','Interpreter','Latex');
drawnow
end

