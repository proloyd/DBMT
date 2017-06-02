function [ final_Confidence_bounds ] = confidence_interval_log_DBMT( n, num_seq )
%[ final_Confidence_bounds ] = confidence_interval_log_DBMT(n, num_seq) 
%plots 95% confidence interval for log_DBMT estimate at nth  
%window assuming num_seq number of tapers were used in that  
%estimate. The multiplier can be changed to obtain any other 
%confidence level.
%
% Inputs:
%   n = window number
%   num_seq = # of tapers used
%
% Outputs:
%   final_Confidence_bounds(:,1) = upper confidence bound
%   final_Confidence_bounds(:,2) = lower confidence bound


multiplier = 2.326; % From normal table
Fs = 110; % in Hz  , Sampling Frequency (Change Accordingly)
W = 6;  % s , Window length (Change Accordingly)
W = W*Fs;
final_est = 0;
final_Confidence_bounds = 0;
pathname = fileparts('./log_DBMT/file');
for j = 1:num_seq
    file_name = sprintf('taper%d.mat',j);
    matfile = fullfile(pathname, file_name);
    load(matfile);
    x_n = x(:,n);
    variance_n = P(:,:,n);
    Confidence_bounds = [x_n + multiplier*sqrt(diag(variance_n)),x_n - multiplier*sqrt(diag(variance_n))];
    % figure, plot(x_n);
    % hold on;
    % plot(Confidence_bounds(:,1),'g');
    % plot(Confidence_bounds(:,2),'r');
    
    final_est = final_est + exp(x_n);
    final_Confidence_bounds = final_Confidence_bounds + exp(Confidence_bounds);

end

final_est = final_est/3;
final_Confidence_bounds = final_Confidence_bounds/3;

%%
plot((0:120)/6,10*log10(final_est(1:121)/Fs));
hold on
plot((0:120)/6,10*log10(final_Confidence_bounds(1:121,1)/Fs),'g');
plot((0:120)/6,10*log10(final_Confidence_bounds(1:121,2)/Fs),'r');
xlim([0 20]);
ylim([-80 20]);
grid on
xlabel('Frequency(Hz)','Interpreter','Latex');
ylabel('PSD(in dB)','Interpreter','Latex');
legend('Est','Upper Bound','Lower Bound','Location','Northwestoutside');
end

