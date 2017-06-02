function [ o_mt_est ] = TSpectrogram( y,N,overlap,time_halfbandwidth,seq_num,Fs )
%[ o_mt_est ] = TSpectrogram( y,N,overlap,time_halfbandwidth,seq_num,Fs )computes 
% single tapered spectrogram. 
%
% Outputs:
%   o_mt_est = conventional multitaper estimate
%
% Inputs:
%   y = Data
%   N = window length in seconds
%   overlap \in [0,1]
%   time_halfbandwidth = N*B
%   num_seq = taper to be used
%   Fs = Sampling Frequency


N = floor(Fs*N);        
K = N*(1-overlap);        
W = floor(length(y)/K - N/K);

seq_length = N;

[dps_seq,lambda] = dpss(seq_length,time_halfbandwidth,seq_num);

o_mt_est = zeros(N,W+1);

for k = 0:1:W
    d = y(1+k*K:k*K+N); 
    o_mt_est(:,k+1) = abs(fft(d.*dps_seq(:,seq_num))).^2;
end

end

