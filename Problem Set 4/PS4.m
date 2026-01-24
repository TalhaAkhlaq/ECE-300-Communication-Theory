%% Talha Akhlaq
% ECE300 Communication Theory
% Problem Set 4
clc; clear; close all;

%% P3

fs = 100;
N = 2000;
t = (0:N-1)/fs;
fc = 10;
fm = 2;
Ac = 2;
ka = 2;
Am_vec = [0.40 0.45 0.60];
idx = 1:400;
for k = 1:length(Am_vec)
    Am = Am_vec(k);
    m = Am*sin(2*pi*fm*t);
    env = Ac*(1 + ka*m);
    x = env.*cos(2*pi*fc*t);

    figure;
    plot(t(idx), x(idx)); hold on;
    plot(t(idx),  env(idx), '--');
    plot(t(idx), -env(idx), '--');
    hold off;
    xlabel('t (s)');
    ylabel('x_{AM}(t)');
    title(sprintf('A_m = %.2f', Am));
end

%% P4

N = 1e6;
rho = sqrt(6*log(10));
sigma2_list = [1 10^(-0.1) 10^(0.1)];

for k = 1:length(sigma2_list)
    s2 = sigma2_list(k);
    s = sqrt(s2);
    nI = s*randn(1,N);
    nQ = s*randn(1,N);
    R  = sqrt(nI.^2 + nQ.^2);

    figure;
    histogram(R, 100, 'Normalization', 'pdf'); hold on;
    yl = ylim;
    plot([rho rho], yl, 'k', 'LineWidth', 2);
    hold off;
    xlabel('r');
    ylabel('pdf estimate');
    if k == 1
        title('\sigma^2 = 1 (0 dB)');
    elseif k == 2
        title('\sigma^2 = 10^{-0.1} (-1 dB)');
    else
        title('\sigma^2 = 10^{0.1} (+1 dB)');
    end
end
