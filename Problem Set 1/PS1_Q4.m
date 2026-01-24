%% Talha Akhlaq
% ECE300 Communication Theory
% Problem Set 1: Communications Overview
% Question 4
clc; clear; close all;

%% Part (b):
S_psk8  = psk_constellation(8);         % 2x8; s_m = A*exp(j*2π(m−1)/8); A = 1/|exp(j*2π/8)−1| => d_min=1
S_qam16 = qam16_constellation();        % 2x16; I,Q ∈ {±1,±3}; scale 1/2 → nearest spacing=1 => d_min=1
S_orth4 = orthogonal_constellation(4);  % 4x4; sqrt(E0)*e_i; d_min=sqrt(2E0); set E0=1/2 => d_min=1

%% Part (a):
[name1, met1] = metrics_named('8-PSK',  S_psk8);
[name2, met2] = metrics_named('16-QAM', S_qam16);
[name3, met3] = metrics_named('4-orth', S_orth4);

%% Part (c):
Names = string({name1; name2; name3});
Mvals = [met1.M; met2.M; met3.M];
Nvals = [met1.N; met2.N; met3.N];
dmin  = [met1.dmin; met2.dmin; met3.dmin];
Es    = [met1.Es; met2.Es; met3.Es];
Eb    = [met1.Eb; met2.Eb; met3.Eb];
EbdB  = [met1.Eb_dB; met2.Eb_dB; met3.Eb_dB];
eta   = [met1.eta; met2.eta; met3.eta];

T = table(Names, Mvals, Nvals, dmin, Es, Eb, EbdB, eta, ...
    'VariableNames', {'Constellation','M','N','d_min','E_s','E_b', ...
                      'E_b_dB_per_bit','eta_bits_per_dimension'});
disp(T);

%% Part (d):
% Most spectrally efficient
[~, idx_eta_max] = max(eta);
fprintf('\nMost spectrally efficient: %s (eta = %.4f bits/dimension)\n', ...
    Names(idx_eta_max), eta(idx_eta_max));

%% Part (e):
% Most power efficient (smallest E_b at fixed d_min)
[~, idx_power_min] = min(Eb);
fprintf('Most power efficient: %s (E_b = %.6f, %.4f dB/bit)\n', ...
    Names(idx_power_min), Eb(idx_power_min), EbdB(idx_power_min));

% dB/bit difference between most and least spectrally efficient
[~, idx_eta_min] = min(eta);
delta_dB = EbdB(idx_eta_max) - EbdB(idx_eta_min);
fprintf('Difference in dB/bit between most and least spectrally efficient: %.4f dB\n', delta_dB);

%% Functions

% PSK constellation
function S = psk_constellation(M)
A = 1/abs(1 - exp(1j*2*pi/M));                 % d_min = A*exp(j*2π/M)−1| = 1
m = 0:M-1;                                     % m = 0,...,M-1
z = A*exp(1j*2*pi*m/M);                        % Complex points s_m = A*exp(j*2πm/M), m=0…M−1
S = [real(z); imag(z)];                        % 2×M; columns = [Re(s_m); Im(s_m)]
assert(abs(min_pairwise_dist(S) - 1) < 1e-12); % Check d_min = 1
end


% 16-QAM constellation (scaled so d_min = 1)
function S = qam16_constellation()
levels = [-3 -1 1 3];                % Amplitude levels (I,Q ∈ {±1, ±3})
[I,Q] = meshgrid(levels, levels);    % Grid of 16 points (I/Q plane)
Sraw = [I(:)'; Q(:)'];               % 2×16 raw constellation matrix
s = 0.5;                             % Scale factor: raw spacing = 2 => scaled so d_min = 1
S = s*Sraw;                          % Final constellation (2×16, spacing normalized)
assert(abs(min_pairwise_dist(S) - 1) < 1e-12);
end

% M-ary orthogonal constellation in R^M, scaled so d_min = 1
function S = orthogonal_constellation(M)
E0  = 1/2;                   % d_min^2 = 2*E0  =>  E0 = 1/2  (with d_min = 1)
amp = sqrt(E0);              % s_i = sqrt(E0)*e_i => amplitude = sqrt(E0)
S   = amp * eye(M);          % constellation matrix: columns are sqrt(E0)*e_i (orthogonal basis vectors in R^M)
assert(abs(min_pairwise_dist(S) - 1) < 1e-10);
end

% Constellation metrics from S ∈ R^{N×M} (columns = symbol vectors)
function [name, met] = metrics_named(name, S)
[N, M] = size(S);
Es = mean(sum(S.^2,1));            % Energy per symbol
Eb = Es / log2(M);                 % Energy per bit 
met = struct();                    % Struct to store constellation metrics
met.N = N; met.M = M; met.Es = Es; met.Eb = Eb;
met.Eb_dB = 10*log10(Eb);          % Eb in dB (true γb would be Eb/N0, but omitted N0 since constellations compared at same d_min)
met.eta = log2(M)/N;               % Spectral efficiency η ≈ (log2 M)/N (bits per real dimension)
met.dmin = min_pairwise_dist(S);
end

% Minimum distance
function d = min_pairwise_dist(S)
[~, M] = size(S);                      % Columns of S = number of symbols M
d = inf;                               % Start with d = infinity (so any distance will be smaller)
for i = 1:M-1
    D = vecnorm(S(:,i) - S(:,i+1:M));  % Distances from symbol i to all later symbols
    d = min(d, min(D));                % Update d with the smallest distance found so far
end
end
