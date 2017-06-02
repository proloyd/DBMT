function [ u_b,l_b ] = confidence_interval_CMT( n, num_seq )
%[ u_b,l_b ] = confidence_interval_CMT( n, num_seq ) computes 
%95% confidence interval of conventional multitaper estimates
%of a single window.
%
% Inputs:
%   n = window number
%   num_seq = # of tapers used
%
% Outputs:
%   u_b = upper confidence bound
%   l_b = lower confidence bound

degrees_of_freedom = 2*num_seq;
u_mul = 1.237;
l_mul = 12.592;
Fs = 110;     % Sampling Frequency  (Change accordingly)
N = 6*Fs;     % 6 sec window size   (Change accordingly)

pathname = fileparts('./CMT/file');
matfile = fullfile(pathname, 'mtm_est.mat');
load(matfile);

S = c_mt_est(:,n);
u_b = (degrees_of_freedom-1)*S/u_mul;
l_b = (degrees_of_freedom-1)*S/l_mul;

plot(Fs*(0:(21*N/Fs)-1)/N,10*log10(S(1:21*N/Fs)/Fs));
hold on;
plot(Fs*(0:(21*N/Fs)-1)/N,10*log10(u_b(1:21*N/Fs)/Fs),'g');
plot(Fs*(0:(21*N/Fs)-1)/N,10*log10(l_b(1:21*N/Fs)/Fs),'r');
xlim([0 20]);
ylim([-80 20]);
grid on
xlabel('Frequency(Hz)','Interpreter','Latex');
ylabel('PSD(in dB)','Interpreter','Latex');
legend('Est','Upper Bound','Lower Bound','Location','Northwestoutside');
end

