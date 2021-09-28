% FFT formulation for aortic thoracic pressure
clear all;
clc;

% from mmgh to Pascal
CF=101325/760;
% Period
T = 60/72;
% Create vector of time
time = T*[0:0.03125:1];
% Create points approximating to [Kim,2010] findings for aortic back pressure
pressure = CF*[78 77 76.5 77 90 115 120.5 120 119.5 119 116 111 103 96 93 93 92 ...
91 90 89 88 87 86 85.25 84.5 83.75 83 82.25 81.5 80.75 80 79.25 78.5];
% Plot these points


plot(time,pressure,'ro','Linewidth',2);
xlim([0 T]);
xlabel('time (seconds)');
ylabel('pressure (Pa)');
%to remove exponential in Ylabel
set(gca, 'YTickLabel', num2str(get(gca, 'YTick').'));
grid on;
% Calculate a fast Fourier transform
% Define parameters
d = fft(pressure);
m = length(pressure);
M = floor((m+1)/2);
% Calculate coefficients
a0 = d(1)/m;
an = 2*real(d(2:M))/m;
alast = d(M+1)/m;
bn = -2*imag(d(2:M))/m;

hold on;
% Linear trig combinations at varying frequencies
t = 0:0.01:2;
n = 1:length(an);
p = a0 + an*cos(2*pi*n'*t/T) ...
+ bn*sin(2*pi*n'*t/T) ...
+ alast*cos(2*pi*(length(an)+1)*t/T);
% plot Fourier approximation of smooth curve between plotted points
plot(t,p,'Linewidth',2);
% Generate legend and title
legend('Data','DFT Interpolant:');
title({'Pressure at descending aorta and its branches for one period'});
% Print out formulation
