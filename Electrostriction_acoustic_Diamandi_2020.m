%% TAMRs in optical fibers, excitation and detection
% Pump: fundamental mode. Probe: fundamental mode. Only R0m
% Luis Alberto Sánchez Domínguez, Laboratory of Fiber Optics
% Ref: Diamandi et al (2017). Diamandi et al (2020). 
%%
close all, clear all

% Fiber parameters
rho0 = 2.203*1e3; % Density of silica [kg/m^3]
c = 299792458; % [m/s]
a = 62.5e-6; % Fiber radius [m]
VL = 5949.57; % [m/s] Longitudinal acoustic velocity for fused silica
d1 = 8.0e-6; % Mode field diameter SMF28@1064 [m] 
w1 = d1/(2*sqrt(2)); % Gaussian beam radius [m] 
nr = 1.45; % Refractive index of silica (at 1064 nm)
P11 = 0.121;
P12 = 0.27;
a1 = -nr^4*(P11-P12);
a2 = -nr^4*P12;

% Pump laser parameters
frep = 19.55e3; % Repetition frequency [Hz]
T0 = 700e-12; % Pulse duration [s]
Pavg = 60e-3; % Average power in fiber [W]
Pp = Pavg/(T0*frep); % Peak power [W]

% Simulation parameters
vWindow = 20e9; % Frequency span [Hz]
n = 2^21;
dv = vWindow/n; % Frequency domain resolution [Hz]
v = (-vWindow/2 : dv : vWindow/2 - dv); % Angular frequency [rad/s]
t = (-n/2 : 1 : (n/2-1))/vWindow;
dt = 1/vWindow;
tWindow = n/vWindow;
%% Temporal properties Pump pulse
Ppump = Pp*sech((t-3e-9)/(T0/1.76)).^2; % Instantaneous optical power [W]
plot(t*1e9,Ppump/max(Ppump))
xlabel('Time [ns]')
ylabel('Optical power [a.u.]')
xlim([0 6])
%printplot = gcf;
%exportgraphics(printplot,'PumpPulse.png','Resolution',150)

Pspec = fftshift(dt*fft(ifftshift(Ppump))); % Safe to always use this because of a bug
Pspeclog = 10*log10((abs(Pspec).^2/max(abs(Pspec).^2)));
plot(v*1e-6,Pspeclog)
xlabel('Frequency [MHz]')
ylabel('Power spectral density [dB]')
xlim([0 1000])
%printplot = gcf;
%exportgraphics(printplot,'PumpSpectrum.png','Resolution',150)
clear Pp Ap Pspeclog

%% Integral Qes
E0 = @(r) 1/(sqrt(pi)*w1)*exp(-r.^2/(2*w1^2)); % Gaussian transverse profile
dE0 = @(r) -r/(sqrt(pi)*w1^3).*exp(-r.^2/(2*w1^2)); % dE0/dr

% figure(2);
% subplot(1,2,1)
% fplot(E0, [-a a])
% hold on
% % fplot(E0Bel, [-a a])
% xlabel('r [m]')
% ylabel('E_0 [a.u.]')
% subplot(1,2,2)
% fplot(N2t_E02,[-2*d 2*d])
% xlabel('r [m]')
% ylabel('\nabla_t {E_0}^2 [a.u.]')

xim = [1.98761995267414,5.37108271665326,8.56060387138268,11.7236152202587,14.8774238490686,18.0269279196606,21.1740665294565,24.3197643164307,27.4645189461084,30.6086223121546,33.7522570164956,36.8955431804631,40.0385630795169,43.1813749680301,46.3240212536744] % 15 first solutions of acoustic modes alfa = 0.6305
freqm = VL*(xim.')/(2*pi*a); % [Hz]

% Gamma_m from bandwidths
gamma_m = [22511.34598 23887.88044 31537.56255 42016.20024 54277.07599 67272.33851 82498.2069 95881.9198 114655.17249 131059.01164 149585.25835 169512.26718 181925.42015 206836.38613 220460.87479];

for m = 1:length(xim)
    um_den = @(r) r.*besselj(1,xim(m)*r/a).^2;
    um_integral = sqrt(2*pi*integral(um_den,0,a));
    um = @(r) besselj(1,xim(m)*r/a)/um_integral;
    fun = @(r) E0(r).*dE0(r).*um(r).*r;
    Qes(m) = -(a1+4*a2)*2*pi*integral(fun,0,a);
end
Qes = Qes'; gamma_m = gamma_m';

A_v = (1./(4*nr*c*rho0*gamma_m.*freqm)).*Qes.*Pspec./(1i - 2*(v - freqm)./gamma_m); % - sign: forward propagation?
%clear Pspec
% figure(3)
% plot(v*1e-6,10*log10(abs(A_v(5,:)).^2/max(abs(A_v(:)).^2)));
% %plot(v*1e-6,(abs(h_w)))
% %xlim([0 1000])
% %ylim([-50 0])
% xlabel('Frequency [MHz]')
% ylabel('Power spectral density (h_m) [dB]')

%% Integral Qpe
for m = 1:length(xim)
    um_den = @(r) r.*besselj(1,xim(m)*r/a).^2;
    um_integral = sqrt(2*pi*integral(um_den,0,a));
    um = @(r) besselj(1,xim(m)*r/a)/um_integral;
    dum = @(r) (xim(m)/2*a).*(besselj(0,xim(m)*r/a) - besselj(2,xim(m)*r/a))/um_integral;
    fun = @(r) (dum(r) + um(r)./r).*abs(E0(r)).^2.*r;
    Qpe(m) = (a1/2 + a2)*2*pi*integral(fun,0,a);
end
Qpe = Qpe';
%% Effective index variation
dn_v = (1./(4*nr^2*c*rho0*gamma_m.*freqm)).*Qes.*Qpe.*Pspec./(1i - 2*(v - freqm)./gamma_m); %  - sign num?
dn_t = (2*pi/sqrt(n))*ifft(ifftshift(dn_v,2),[],2);
sumadn=sum(dn_t,1);
clear dneff
figure(5)
plot((t+abs(min(t)))*1e6,real(sumadn),'black')
xlim([0 3])

% Dont apply a hanning, hamming etc window to a damped signal
Pspecdneff = fftshift(dt*fft(ifftshift(sumadn)));
%Pspecdneff = fftshift(dt*fft(ifftshift(sumadn(1000:end)))); % Quitamos el pico inicial

figure(6)
plot(v*1e-6,10*log10((abs(Pspecdneff).^2)/max(abs(Pspecdneff).^2)),'black')
%plot(10*log10((abs(Pspecdneff).^2)/max(abs(Pspecdneff).^2)),'black') %Quitamos el pico inicial
xlim([0 1200])
xlabel('Frequency (MHz)','FontSize',16)
ylabel('Relative intensity (dB)','FontSize',16)
set(gca,'FontSize',14);
