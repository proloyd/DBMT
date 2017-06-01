function [ final_Confidence_bounds ] = confidence_interval_DBMT( n, num_seq )
%confidence_interval_DBMT(n, num_seq) plots 95% confidence interval for
%DBMT estimate at nth window assuming num_seq number of tapers were used 
%in that estimate. The multiplier can be changed to obtain any other 
%confidence level.

multiplier = 2.326; % From normal table
Fs = 110; % in Hz
W = 6;  % s 
W = W*Fs;
final_est = 0;
final_Confidence_bounds = 0;
pathname = fileparts('./DBMT/file');
for j = 1:num_seq
    file_name = sprintf('taper%d.mat',j);
    matfile = fullfile(pathname, file_name);
    load(matfile);
    x_n = x(:,n);
    variance_n = P(:,:,n);
    Confidence_bounds = [x_n + multiplier*sqrt(diag(variance_n)),x_n - multiplier*sqrt(diag(variance_n))];
    % figure, plot(x);
    % hold on;
    % plot(Confidence_bounds(:,1),'g');
    % plot(Confidence_bounds(:,2),'r');

    Confidence_bounds1 = zeros(size(Confidence_bounds));
    for i = 1:length(x_n)
        if Confidence_bounds(i,1)*Confidence_bounds(i,2) < 0
            Confidence_bounds1(i,:) = [max(Confidence_bounds(i,:).^2),10^-15];
        else
            Confidence_bounds1(i,:) = sort((Confidence_bounds(i,:)).^2,'descend');
        end
    end
    % figure, plot(x.^2);
    % hold on;
    % plot(Confidence_bounds1(:,1),'g');
    % plot(Confidence_bounds1(:,2),'r');

    est = (x_n(1:length(x_n)/2).^2 + flipud(x_n(length(x_n)/2+1:length(x_n))).^2)/4;
    Confidence_bounds2 = (Confidence_bounds1(1:length(x_n)/2,:) + flipud(Confidence_bounds1(1+length(x_n)/2:end,:)))/4;
    % figure, plot(est);
    % hold on
    % plot(Confidence_bounds2(:,1),'g');
    % plot(Confidence_bounds2(:,2),'r');
    
    final_est = final_est + est;
    final_Confidence_bounds = final_Confidence_bounds + Confidence_bounds2;
end
final_est = final_est/3;
final_Confidence_bounds = final_Confidence_bounds/3;

%%
plot((1:120)/6,10*log10(W^2*final_est((1:120))/Fs));
hold on
plot((1:120)/6,10*log10(W^2*final_Confidence_bounds(1:120,1)/Fs),'g');
plot((1:120)/6,10*log10(W^2*final_Confidence_bounds(1:120,2)/Fs),'r');
xlim([0 20]);
ylim([-80 20]);
grid on
xlabel('Frequency(Hz)','Interpreter','Latex');
ylabel('PSD(in dB)','Interpreter','Latex');
legend('Est','Upper Bound','Lower Bound','Location','Northwestoutside');

end

