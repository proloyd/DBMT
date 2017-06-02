function [ c_mt_est ] = MTSpectrogram( y,W,overlap,time_halfbandwidth,num_seq,Fs )
%MTSpectogram( y,N,overlap,time_halfbandwidth,num_seq,Fs ) computes overlapped multitaper
% Spectrogram 

%   Outputs:
%   c_mt_est = conventional multitaper estimate
%
%   Inputs:
%   y = Data
%   W = window length in seconds
%   overlap \in [0,1]
%   time_halfbandwidth = N*B
%   num_seq = # of tapers to be used
%   Fs = Sampling Frequency

W = Fs*W;        
K = W*(1-overlap);        
N = floor(length(y)/K - W/K);

seq_length = W;
[dps_seq,lambda] = dpss(seq_length,time_halfbandwidth,num_seq);

c_mt_est = zeros(W,N+1);

for k = 0:1:N
    d = y(1+k*K:k*K+W); 
    for j = 1:num_seq
        c_mt_est(:,k+1) = c_mt_est(:,k+1)+ abs(fft(d.*dps_seq(:,j))).^2;
    end
end
c_mt_est = c_mt_est/num_seq;
%****************************Plot upto 20 Hz*******************************
pcolor(((0:1:N)*K+W/2)/Fs,Fs*(0:(21*W/Fs)-1)/W,10*log10(c_mt_est(1:21*W/Fs,:)/Fs));
shading flat;
colormap('jet');
colorbar;
xlabel('Time(s)','Interpreter','Latex');
ylabel('Frequency(Hz)','Interpreter','Latex');
drawnow;
title('Overlapped Multi-taper Spectrogram','Interpreter','Latex');
pathname = fileparts('./CMT/file');
matfile = fullfile(pathname, 'mtm_est');
save(matfile, 'c_mt_est');
end

