%% Hüseyin Berk Keskin EEE409 Project 2
clc, clear, close all;

% System Parameters
M = 2;                          % modulation order
Tb = 1;                         % bit interval(sec)
ts = 0.01;                      % sampling times of the pulse (for Matlab realization)
Nbits = 10^6;                   % number of transmitted bits
Ac = 10;                        % Carrier amplitude for binary input '1'
fc = 10 / Tb;                   % Carrier frequency
SNR_dB = 0:12;                  % SNR range in dB
Eb = (Ac^2 * Tb) / 2;           % Energy per bit
SNR_lin = 10.^(SNR_dB / 10);    % Convert SNR from dB to linear scale
phase_offset = pi / 11.5;       % Phase offset for Bob

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Binary data and secret key generation
data = randi([0, 1], 1, Nbits);         % Generate random binary data
secret_key = randi([0, 1], 1, Nbits);  % Generate secret key for encryption

% XOR for encryption
encrypted_data = xor(data, secret_key);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Preallocate BER arrays
BOB_BER = zeros(size(SNR_dB));
EVE_BER = zeros(size(SNR_dB));
THE_BER = zeros(size(SNR_dB));
PLS_BER_BOB = zeros(size(SNR_dB));
PLS_BER_EVE = zeros(size(SNR_dB));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Modulate encrypted data using Binary ASK
carr_sig = Ac * cos(2 * pi * fc * (0:ts:Tb - ts)) * sqrt(2 / Tb);
mod_sig = reshape(repmat(encrypted_data, length(carr_sig), 1) .* carr_sig', 1, []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Iterate over each SNR value
for snr_idx = 1:length(SNR_dB)
    SNR = SNR_lin(snr_idx);
    N0 = Eb / SNR;
    noise_std_dev = sqrt(N0 / 2);

    % Add noise to Bob's received signal
    BOB_NOISE = noise_std_dev * randn(1, length(mod_sig));
    BOB_RECEIVED = mod_sig + BOB_NOISE;

    % Add noise and distortions to Eve's received signal
    EVE_rand_phase = 2 * pi * rand(1, length(mod_sig));
    EVE_rand_amplitude = 0.5 + 0.5 * rand(1, length(mod_sig));
    EVE_dist_sig = EVE_rand_amplitude .* mod_sig .* cos(EVE_rand_phase);
    EVE_NOISE = noise_std_dev * randn(1, length(EVE_dist_sig));
    EVE_RECEIVED = EVE_dist_sig + EVE_NOISE;

    % Bob's demodulation
    BOB_demod_data = arrayfun(@(idx) sum(BOB_RECEIVED(idx:idx + length(carr_sig) - 1) .* carr_sig) * ts > 0.5 * Ac, 1:length(carr_sig):length(mod_sig));
    BOB_decrypted_data = xor(BOB_demod_data, secret_key);
    BOB_BER(snr_idx) = mean(BOB_decrypted_data ~= data);

    % Eve's demodulation
    EVE_demod_data = arrayfun(@(idx) sum(EVE_RECEIVED(idx:idx + length(carr_sig) - 1) .* carr_sig) * ts > 0.5 * Ac, 1:length(carr_sig):length(mod_sig));
    EVE_BER(snr_idx) = mean(EVE_demod_data ~= encrypted_data);

    % Calculate theoretical BER
    THE_BER(snr_idx) = qfunc(sqrt(2 * SNR));

    % PLS BER method for Bob and Eve
    TX_data = 2 * data - 1; % Map binary data to +/-1
    modulated_signal = (TX_data + 1) .* exp(-phase_offset);
    noise_variance = 1 / SNR_lin(snr_idx);
    sigma = sqrt(noise_variance);
    noise = sigma * randn(1, Nbits);
    RX_data = modulated_signal + noise;

    RX_Bob = RX_data .* exp(phase_offset) - 1; % Correct Bob's phase
    RX_Eve = RX_data; % No correction for Eve

    detected_Bob = sign(RX_Bob);
    detected_Eve = sign(RX_Eve);

    PLS_BER_BOB(snr_idx) = sum(0.5 * abs(TX_data - detected_Bob)) / Nbits;
    PLS_BER_EVE(snr_idx) = sum(0.5 * abs(TX_data - detected_Eve)) / Nbits;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot results
figure;
semilogy(SNR_dB, BOB_BER, 'b-o', 'LineWidth', 1.5, 'DisplayName', 'BOB');
hold on;
semilogy(SNR_dB, EVE_BER, 'r--s', 'LineWidth', 1.5, 'DisplayName', 'EVE');
semilogy(SNR_dB, THE_BER, 'k-*', 'LineWidth', 1.5, 'DisplayName', 'THEORETICAL');
semilogy(SNR_dB, PLS_BER_BOB, 'g-^', 'LineWidth', 1.5, 'DisplayName', 'PLS BER BOB');
semilogy(SNR_dB, PLS_BER_EVE, 'm-v', 'LineWidth', 1.5, 'DisplayName', 'PLS BER EVE');
grid on;
xlabel('SNR (dB)');
ylabel('BER');
title('BER Performance Comparison: ASK vs PLS BER');
legend('show');

% Print results
fprintf('SNR(dB)  ASK BOB  ASK EVE  BER BOB  BER EVE  THEORETICAL\n');
for snr_idx = 1:length(SNR_dB)
    fprintf('%d\t\t %.5f  %.5f  %.5f  %.5f   %.5f\n', SNR_dB(snr_idx), BOB_BER(snr_idx), EVE_BER(snr_idx), PLS_BER_BOB(snr_idx), PLS_BER_EVE(snr_idx), THE_BER(snr_idx));
end