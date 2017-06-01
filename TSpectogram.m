function [ o_mt_est ] = TSpectogram( y,N,overlap,time_halfbandwidth,seq_num,Fs )
%UNTITLED11 Summary of this function goes here
%   N = window length in seconds
%   overlap \in [0,1]

N = floor(Fs*N);        % 6 sec window
K = N*(1-overlap);        % 2 sec overlap
W = floor(length(y)/K - N/K);
%%
seq_length = N;
% time_halfbandwidth = 3;
% num_seq = floor(2*(3)-2);
[dps_seq,lambda] = dpss(seq_length,time_halfbandwidth,seq_num);
%%
o_mt_est = zeros(N,W+1);

for k = 0:1:W
    d = y(1+k*K:k*K+N); 
%     o_mt_est(:,k+1) = abs(fft(d.*dps_seq(:,seq_num))).^2;
%     o_mt_est(:,k+1) = abs(fft(d)).^2;
    o_mt_est(:,k+1) = abs(fft(d.*dps_seq(:,seq_num))).^2;
%     for j = 1:num_seq
%         o_mt_est(:,k+1) = o_mt_est(:,k+1)+ abs(fft(d.*dps_seq(:,j))).^2;
%     end
end

% plot upto 20 Hz
% figure, 
% subplot(2,1,1); pcolor(((0:1:W)*K+N/2)/Fs,Fs*(0:(20*N/Fs)-1)/N,10*log10(o_mt_est(1:20*N/Fs,:)/Fs));
% shading flat;
% colormap('jet');
% colorbar;
% % caxis([-40 10])
% xlabel('Time(s)','Interpreter','Latex');
% ylabel('Frequency(Hz)','Interpreter','Latex');
% title('MTSpecrogram - 6s Window, 3s overlap','Interpreter','Latex');

end

