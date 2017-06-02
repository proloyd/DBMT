%Dynamic Bayesian Multi Taper Estimation algorithms
%Master script to regenerate the Fig. 3 in Dynamic Bayesian Multitaper  
%Spectral Analysis.

clear all
mkdir CMT
mkdir DBMT
mkdir log_DBMT
%*************************Generation of Toy Example***********************
T = 600; % in s
Fs = 110;
f0 = 0.02; % in Hz
f1 = 11; %in Hz
f2 = 5; % in Hz;
t = 0:(T/(Fs*600-1)):T;
sigma_b = 6*10^-4;
sigma_a = 1*10^-5;

%************************Amplitude Modulated Component*********************
b = conv(conv(conv([1 -0.98*exp(1j*2*(pi)*f1/110)],[1 -0.98*exp(-1j*2*(pi)*f1/110)]),conv([1 -0.98*exp(1j*2*(pi)*f1/110)],[1 -0.98*exp(-1j*2*(pi)*f1/110)])),conv([1 -0.95*exp(1j*2*(pi)*f1/110)],[1 -0.95*exp(-1j*2*(pi)*f1/110)]));
y1 = filter(1,b,sigma_b*randn(size(t)));
y1 = y1(end-66000+1:end).*((cos(2*pi*f0*t)).^8+0.1);

%************************Frequency Modulated Component*********************
f2d = f2-12/25;
a = conv(conv(conv([1 -0.95*exp(1j*(2*pi)*f2d/110)],[1 -0.95*exp(-1j*2*(pi)*f2d/110)]),conv([1 -0.98*exp(1j*2*(pi)*f2d/110)],[1 -0.98*exp(-1j*2*(pi)*f2d/110)])),conv([1 -0.98*exp(1j*2*(pi)*f2d/110)],[1 -0.98*exp(-1j*2*(pi)*f2d/110)]));
c = conv(conv([1 1],[1 1]),conv([1 1],[1 1]));
y2 = filter(c,a,0.1*randn(5*110,1));
for k = 1:25
    f2d = f2d+12/25;
    a = conv(conv(conv([1 -0.95*exp(1j*(2*pi)*f2d/110)],[1 -0.95*exp(-1j*2*(pi)*f2d/110)]),conv([1 -0.98*exp(1j*2*(pi)*f2d/110)],[1 -0.98*exp(-1j*2*(pi)*f2d/110)])),conv([1 -0.98*exp(1j*2*(pi)*f2d/110)],[1 -0.98*exp(-1j*2*(pi)*f2d/110)]));
    y_inter = filter(c,a,sigma_a*randn(26*110,1)*(1.17)^(k-1));
    y2 = [y2;y_inter(221:26*110)];
end
y2 = y2(end-66000+1:end);

%******************************True data y_true****************************
y_true = y1' + y2;
% figure, plot(t,y_true,'k');
% xlabel('time','Interpreter','Latex');
% ylabel('Amplitude','Interpreter','Latex');

%**************************Noisy version of data, y************************
sigma_n = sqrt(var(y_true)/10^(30/10));
y = y_true + sigma_n*randn(size(y_true));

figure, plot(t,y,'k');
title('Noisy AR process');
xlabel('time','Interpreter','Latex');
ylabel('Amplitude','Interpreter','Latex')

%****************************Ground Truth***********************************
sigma_b =(5*10^-4);
[Hb,Freq] = freqz(1,b,330,Fs);
% figure, plot(Freq, 20*log10(abs(Hb)));
var_b = cos(2*pi*f0*t).^8;
PSDb = [];
for n = 1:100
    PSDb = [PSDb, (sigma_b*(mean(var_b((n-1)*660+1:n*660))+0.1).^2)*(Hb)];
end
f2d = f2-12/25;
a = conv(conv(conv([1 -0.95*exp(1j*(2*pi)*f2d/110)],[1 -0.95*exp(-1j*2*(pi)*f2d/110)]),conv([1 -0.98*exp(1j*2*(pi)*f2d/110)],[1 -0.98*exp(-1j*2*(pi)*f2d/110)])),conv([1 -0.98*exp(1j*2*(pi)*f2d/110)],[1 -0.98*exp(-1j*2*(pi)*f2d/110)]));
c = conv(conv([1 1],[1 1]),conv([1 1],[1 1]));
PSDa = [];
for k = 1:25
    f2d = f2d+12/25;
    a = conv(conv(conv([1 -0.95*exp(1j*(2*pi)*f2d/110)],[1 -0.95*exp(-1j*2*(pi)*f2d/110)]),conv([1 -0.98*exp(1j*2*(pi)*f2d/110)],[1 -0.98*exp(-1j*2*(pi)*f2d/110)])),conv([1 -0.98*exp(1j*2*(pi)*f2d/110)],[1 -0.98*exp(-1j*2*(pi)*f2d/110)]));
    [Ha,Freq] = freqz(c,a,330,Fs);
    PSDa = [PSDa, (1.17)^(k-1)*(sigma_a*(Ha)),(1.17)^(k-1)*(sigma_a*(Ha)),(1.17)^(k-1)*(sigma_a*(Ha)),(1.17)^(k-1)*(sigma_a*(Ha))];
end
PSD_y = 10*log10((abs(PSDa+PSDb).^2+sigma_n^2)/(Fs));
PSD_true = 10*log10((abs(PSDa+PSDb).^2)/(Fs));
figure, 
subplot(4,1,1), pcolor(((0:99)'+0.5)*6,Freq(1:300),PSD_y(1:300,1:100));
ylim([0 20]);
shading flat;
colormap('jet');
colorbar;
xlabel('Time(s)','Interpreter','Latex');
ylabel('Frequency(Hz)','Interpreter','Latex');
title('Ground Truth','Interpreter','Latex');
drawnow
%*****************************Estimates************************************
W = 6;          % 6 s window-length
R = Fs*6;       % freq resolution
U = 329;        % # of frequency bins to be estimated
N = floor(length(y)/(Fs*W));
rho = 3;
K = 3;          % # of tapers

%************************Multitaper Spectrogram****************************
overlap = 0.5;  % 50% Overlap
subplot(4,1,2),
mtm_est = MTSpectrogram(y,W,overlap,rho,K,Fs);

%***************************DBMT Spectrogram*******************************
subplot(4,1,3),
DBMT_est = DBMTSpectrogram(y,W,R,U,rho,K,Fs);

%*************************log_DBMT Spectrogram*****************************
subplot(4,1,4),
log_DBMT_est = log_DBMTSpectrogram(y,W,U,rho,K,Fs);

%*****Confidence Intervals for window starting at t = 79*6 = 477 sec******* 
figure, 
subplot(4,1,1),
plot(Freq(1:150),PSD_y(1:150,79))
hold on
plot(Freq(1:150),PSD_true(1:150,79))
xlim([0 20]);
ylim([-80 20]);
xlabel('Frequency(Hz)','Interpreter','Latex');
ylabel('PSD(in dB)','Interpreter','Latex');
legend('Noisy','True','Location','Northwestoutside');
grid on
subplot(4,1,2)
confidence_interval_CMT(157,3);
plot(Freq(1:150),PSD_y(1:150,79),'k-.');    %reference
subplot(4,1,3)
confidence_interval_DBMT(79,3);
plot(Freq(1:150),PSD_true(1:150,79),'k-.'); %reference
subplot(4,1,4)
confidence_interval_log_DBMT(79,3);
plot(Freq(1:150),PSD_y(1:150,79),'k-.');    %reference
