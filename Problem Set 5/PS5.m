%% Talha Akhlaq
% ECE300 Communication Theory
% Problem Set V: Digital Communications
clear; close all; clc;

q = @(x) 0.5 * erfc(x / sqrt(2));

%% Problem 1

% SNR: 4dB
g = 10^(4/10); 

% Coherent Orthogonal (rho = 0)
coh_orth = q(sqrt(g));
fprintf('Problem 1:\n');
fprintf('\tCoherent Orthogonal: %.4f\n', coh_orth);

% Coherent Optimal FSK (rho = -0.217)
coh_opt = q(sqrt((1 - (-0.217)) * g)); 
fprintf('\tCoherent Optimal FSK: %.4f\n', coh_opt);

%% Problem 2 (d)

a = [2.5, 1.5, 1.0, 1.0, 0.5, 0.5];
b = [1,   2,   4,   5,   9,   10 ];

gamma_dB = 4:0.1:14;
g = 10.^(gamma_dB / 10);

Pe_detailed = zeros(size(g));
Pe_simple   = zeros(size(g));

for i = 1:length(g)
    Pe_detailed(i) = sum(a .* q(sqrt(b .* g(i))));
    Pe_simple(i)   = a(1) * q(sqrt(b(1) * g(i)));
end

idx = find(abs(gamma_dB - 10) < 1e-5);
fprintf('Problem 2(d):\n')
fprintf('\t10 dB:\n\t Detailed Pe = %.4e\n\t Simple Pe   = %.4e\n', ...
    Pe_detailed(idx), Pe_simple(idx));

%% Problem 2 (e) 

figure;
semilogy(gamma_dB, Pe_detailed, 'b-', 'LineWidth', 2); hold on;
semilogy(gamma_dB, Pe_simple, 'r--', 'LineWidth', 2);
yline([1e-2 1e-6], 'w:', 'LineWidth', 1.5); 

xlabel('SNR per bit, \gamma_b (dB)');
ylabel('Probability of Symbol Error, P_e (log)');
title('Problem 2(e): Probability of Error (Pe) vs. SNR per Bit (\gamma_b)');
legend('Detailed Union Bound', 'Simplified (d_{min} only)');
ylim([1e-7 1]); 
xlim([4 14]);
grid on;
set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off');

%% Problem 6 (d) 

% Filter Parameters
rolloff = 0.2;      % Beta
span    = 3;        % Spans 3 symbols
sps     = 16;       % Samples per symbol
shape   = 'sqrt';   % Square-root raised cosine

% Transmit Filter 
h = rcosdesign(rolloff, span, sps, shape);

% Matched Filter Output 
g = conv(h, h);

figure;
% Transmit Pulse
subplot(2,1,1);
stem(h, 'filled');
title('Problem 6 (d): Transmit Pulse (\surdRC)');
xlabel('Sample n'); ylabel('Amplitude');
grid on;

% Matched Filter
subplot(2,1,2);
stem(g, 'filled');
title('Problem 6 (d): Matched Filter (RC)');
xlabel('Sample n'); ylabel('Amplitude');
grid on;

%% Problem 6 (e) 

% Reconstruct Signal
rolloff = 0.2; span = 3; sps = 16;
h = rcosdesign(rolloff, span, sps, 'sqrt');
g = conv(h, h);

% Peak and Interference
[peak_val, center_idx] = max(abs(g));
symbol_indices = center_idx : sps : length(g);
symbol_indices = [fliplr(center_idx:-sps:1), symbol_indices(2:end)];
interference_indices = symbol_indices(symbol_indices ~= center_idx);

% SIR
P_sig = peak_val^2;
I_worst = sum(abs(g(interference_indices))); 
P_int = I_worst^2;                           
SIR0 = 10 * log10(P_sig / P_int);

fprintf('Problem 6(e):\n');
fprintf('\tSignal Power: %.4f\n', P_sig);
fprintf('\tInterference Power: %.4e\n', P_int);
fprintf('\tTheoretical SIR0: %.2f dB\n', SIR0);

%% Problem 6 (f) 

SNIR_target_dB = SIR0 - 5;

% dB to linear
SIR_lin  = 10^(SIR0 / 10);
SNIR_lin = 10^(SNIR_target_dB / 10);

% SNR (linear)
SNR_inv = (1 / SNIR_lin) - (1 / SIR_lin);
SNR_lin = 1 / SNR_inv;

% dB
SNR_required_dB = 10 * log10(SNR_lin);

fprintf('Problem 6(f):\n');
fprintf('\tTarget SNIR: %.2f dB\n', SNIR_target_dB);
fprintf('\tRequired SNR: %.2f dB\n', SNR_required_dB);

%% Problem 7 (a) 

% QPSK Constellation (Gray Coded)
map = [1+1j, -1+1j, -1-1j, 1-1j]; 

% Signal
num_bits = 1e5;
symbol_indices = randi([0 3], 1, num_bits/2); 
syms = map(symbol_indices + 1);              

% Upsample & Filter
syms_up = upsample(syms, sps);
tx_signal = conv(syms_up, h);

% Plot Envelope
figure;
t_plot = 1:1000; 
plot(t_plot, abs(tx_signal(t_plot)), 'LineWidth', 1.5);
title('Problem 7(a): Envelope of Transmitted Waveform');
xlabel('Sample Index'); ylabel('Magnitude'); grid on;

% Average SIR (Transmitted)
[peak_h, center_h] = max(abs(h));
P_sig_tx = 2 * peak_h^2; 

% Interference Power
int_idx = [center_h-sps : -sps : 1,  center_h+sps : sps : length(h)];
P_int_tx = 2 * sum(abs(h(int_idx)).^2);

SIR_tx = 10 * log10(P_sig_tx / P_int_tx);

fprintf('Problem 7(a):\n');
fprintf('\tTransmitted SIR: %.2f dB\n', SIR_tx);

%% Problem 7 (b) 

% Noise Variance 
g = conv(h, h);
[peak_g, center_g] = max(abs(g));

% Signal Power (QPSK Power = 2)
P_sig_out = 2 * peak_g^2;

% Interference Power 
int_idx_g = [center_g-sps : -sps : 1,  center_g+sps : sps : length(g)];
P_int_out = 2 * sum(abs(g(int_idx_g)).^2);

% Noise Power 
if ~exist('SNR_required_dB', 'var'), SNR_required_dB = 9.12; end
P_noise_out_target = P_sig_out / 10^(SNR_required_dB / 10);

% Noise Variance 
filter_energy = sum(abs(h).^2);
sigma2 = P_noise_out_target / filter_energy;

% Noise  
noise = sqrt(sigma2/2) * (randn(size(tx_signal)) + 1j*randn(size(tx_signal)));
rx_noisy = tx_signal + noise;

% Matched Filter 
rx_filtered = conv(rx_noisy, h);

% Plot Envelope 
figure;
t_plot = 1:1000;
plot(t_plot, abs(rx_filtered(t_plot)), 'LineWidth', 1.5);
title('Problem 7(b): Envelope at Matched Filter');
xlabel('Sample Index'); ylabel('Magnitude'); grid on;

% Average SNIR 
Total_Disturbance = P_int_out + P_noise_out_target;
SNIR_out = 10 * log10(P_sig_out / Total_Disturbance);

fprintf('Problem 7(b):\n');
fprintf('\tTarget Output SNR: %.2f dB\n', SNR_required_dB);
fprintf('\tOutput Signal Power: %.4f\n', P_sig_out);
fprintf('\tOutput ISI Power:    %.4f\n', P_int_out);
fprintf('\tOutput Noise Power:  %.4f\n', P_noise_out_target);
fprintf('\tOverall SNIR:        %.2f dB\n', SNIR_out);

%% Problem 7 (c) 
gray_map = [0 0; 0 1; 1 1; 1 0]; 
bits_matrix = gray_map(symbol_indices + 1, :); 
bits = reshape(bits_matrix.', 1, []); 

% Subsample
total_delay = length(h) - 1; 
start_idx = total_delay + 1; 

% Sampling Indices 
num_syms = length(bits) / 2;
sample_indices = start_idx : sps : (start_idx + (num_syms-1)*sps);
max_idx = min(length(sample_indices), length(rx_filtered));
sample_indices = sample_indices(1:max_idx);

% Soft Symbols
rx_soft = rx_filtered(sample_indices);

% Decode Bits
rx_bits = zeros(size(bits));
for i = 1:length(rx_soft)
    r_real = real(rx_soft(i));
    r_imag = imag(rx_soft(i));
    
    if r_real > 0       
        if r_imag > 0      
            b1=0; b2=0;
        else                
            b1=1; b2=0;
        end
    else 
        if r_imag > 0       
            b1=0; b2=1;
        else                
            b1=1; b2=1;
        end
    end
    
    % Store bits
    rx_bits(2*i-1) = b1;
    rx_bits(2*i)   = b2;
end

% Count Errors
bit_errors = sum(bits ~= rx_bits);
BER = bit_errors / length(bits);

% Symbol Errors 
tx_pairs = reshape(bits, 2, []).';
rx_pairs = reshape(rx_bits, 2, []).';
sym_errors = sum(any(tx_pairs ~= rx_pairs, 2));
SER = sym_errors / num_syms;

fprintf('Problem 7(c):\n');
fprintf('\tTotal Bit Errors:    %d\n', bit_errors);
fprintf('\tBit Error Rate:      %.2e\n', BER);
fprintf('\tTotal Symbol Errors: %d\n', sym_errors);
fprintf('\tSymbol Error Rate:   %.2e\n', SER);
