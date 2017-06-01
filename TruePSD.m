T = 600; % s
f0 = 0.02; % Hz
f1 = 5; %Hz;
f2 = 11; %Hz
Fs = 110; %Hz
t = 0:(T/(Fs*600-1)):T;
%%
b = conv(conv(conv([1 -0.99*exp(1j*2*(pi)/10)],[1 -0.99*exp(-1j*2*(pi)/10)]),conv([1 -0.99*exp(1j*2*(pi)/10)],[1 -0.99*exp(-1j*2*(pi)/10)])),conv([1 -0.999*exp(1j*2*(pi)/10)],[1 -0.999*exp(-1j*2*(pi)/10)]));
sigma_b =(5*10^-4);
[Hb,Freq] = freqz(1,b,330,Fs);
% figure, plot(Freq, 20*log10(abs(Hb)));
var_b = cos(2*pi*f0*t).^8;
PSDb = [];
for n = 1:100
    PSDb = [PSDb, (sigma_b*(mean(var_b((n-1)*660+1:n*660))+0.1).^2)*(Hb)];
end
%%
f2d = f1-12/25;
a = conv(conv(conv([1 -0.99*exp(1j*(2*pi)*f2d/110)],[1 -0.99*exp(-1j*2*(pi)*f2d/110)]),conv([1 -0.99*exp(1j*2*(pi)*f2d/110)],[1 -0.99*exp(-1j*2*(pi)*f2d/110)])),conv([1 -0.99*exp(1j*2*(pi)*f2d/110)],[1 -0.99*exp(-1j*2*(pi)*f2d/110)]));
c = conv(conv([1 1],[1 1]),conv([1 1],[1 1]));
sigma_a = 1*10^-5;
PSDa = [];
for k = 1:25
    f2d = f2d+12/25;
    a = conv(conv(conv([1 -0.99*exp(1j*(2*pi)*f2d/110)],[1 -0.99*exp(-1j*2*(pi)*f2d/110)]),conv([1 -0.99*exp(1j*2*(pi)*f2d/110)],[1 -0.99*exp(-1j*2*(pi)*f2d/110)])),conv([1 -0.99*exp(1j*2*(pi)*f2d/110)],[1 -0.99*exp(-1j*2*(pi)*f2d/110)]));
    [Ha,Freq] = freqz(c,a,330,Fs);
    PSDa = [PSDa, (1.17)^(k-1)*(sigma_a*(Ha)),(1.17)^(k-1)*(sigma_a*(Ha)),(1.17)^(k-1)*(sigma_a*(Ha)),(1.17)^(k-1)*(sigma_a*(Ha))];
end

%% 
% fig = 20*log10(abs(PSDb));
% figure, pcolor((1:100)',Freq(1:300),fig(1:300,1:100));
% shading flat;
% colormap('jet');
% colorbar;
% 
% fig = 20*log10(abs(PSDa));
% figure, pcolor((1:100)',Freq(1:300),fig(1:300,1:100));
% shading flat;
% colormap('jet');
% colorbar;

fig = 20*log10(abs(PSDa+PSDb)/sqrt(Fs));
figure, pcolor(((0:99)'+0.5)*6,Freq(1:300),fig(1:300,1:100));
shading flat;
colormap('jet');
colorbar;
%%
figure, plot(Freq(1:150),fig(1:150,15))
xlim([0 20])