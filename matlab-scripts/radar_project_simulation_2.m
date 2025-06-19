clear
clc

%% Parameters

% constants & system parameters
e_0 = 8.8541e-12;           % [F/m]
mu_0 = 1.2566e-6;           % [N/(A^2)]
e_r = 1.00054;
mu_r = 1;
sigma = 1e-9;               % conductivity [S/m]
epsilon = e_r * e_0;
mu = mu_r * mu_0;
c = 1/sqrt(mu * epsilon);

% freq parameters
BW = 4e9;                   % [Hz]
fs = 2 * BW;                % sample rate

f_start = 1e9;        % [Hz]
f_end = 3e9;          % [Hz]

% inputs
N = 2048;               % window size [samples]

% time & freq vectors
t = 0:1/fs:1/fs*(N-1);
f = linspace(0, BW, N);

% Find indices where freq is within the specified range
cut_indices = (f >= f_start) & (f <= f_end);

% Calculate effective bandwidth and number of used frequency points
B = f_end - f_start;   % [Hz]
n = sum(cut_indices);  % number of samples in band

% Frequency step in used region
delta_f = B / n;

% Range resolution and max range
r_min = c / (2 * B);            % [m]
r_max = c / (2 * delta_f);       % [m]

%% amplitude & phase
% Calculate amplitude and phase
[amp, phase] = generate_amp_n_phase(f, sigma, mu, epsilon, r_max);

% Plot the original and filtered amplitude vs frequency
figure(1);
sgtitle(sprintf('Wave propagation (air), z = %.2f(m)', r_max))
subplot(2, 1, 1);
plot(f(cut_indices)*1e-9, 20*log(amp(cut_indices)));
xlabel('Frequency (GHz)');
ylabel('Amplitude (dB)');
title('Amplitude vs Frequency');

subplot(2, 1, 2);
plot(f(cut_indices)*1e-9, phase(cut_indices));
xlabel('Frequency (GHz)');
ylabel('Phase (radians)');
title('Phase vs Frequency');


%% Create a single figure with 3 subplots


% inputs
r_values = [100];   % target distances [m]
% r_values = linspace(r_min, r_max, 4);   % target distances [m]
% r_values = r_values(2:end-1);
figure(2);

for i = 1:length(r_values)
    r = r_values(i);
    
    % Calculate amplitude and phase
    [amp, phase] = generate_amp_n_phase(f, sigma, mu, epsilon, r);
    
    % Compute frequency response and impulse response
    freq_response = amp .* exp(1j * phase);
    freq_response(~cut_indices) = 0;
    imp_response = ifft(freq_response, 'symmetric');
    mag = abs(imp_response);
    
    % Find index of the peak's center
    window_size = 11;
    filtered_mag = movmean(mag, window_size);
    [~, idx_max] = max(filtered_mag);
    detection = zeros(size(mag));
    detection(idx_max) = 1;

    % Plot
    x_axis = t.* c;
    subplot(length(r_values), 1, i);
    plot(x_axis, mag, 'b');  % Impulse response
    hold on;

    % Mark detections
    % plot(x_axis, filtered_mag, 'r');
    plot(x_axis(detection == 1) , mag(detection == 1), 'ro', 'MarkerFaceColor', 'r');
    hold off;

    xlabel('Distance from target (m)');
    ylabel('Amplitude');
    title(sprintf('Impulse Response (r = %.3f m)', r));
    
end


%% Multiple targets simulation


% inputs
r_values = [10, 5, 50, 40, 70];   % target distances [m]
mags = zeros(size(t));

figure(3);
for i = 1:numel(r_values)
    % Calculate amplitude and phase
    [amp, phase] = generate_amp_n_phase(f, sigma, mu, epsilon, r_values(i));
    
    % Compute frequency response and impulse response
    freq_response = amp .* exp(1j * phase);
    freq_response(~cut_indices) = 0;
    imp_response = ifft(freq_response, 'symmetric');
    mag = abs(imp_response);
    
    mags = mags + mag;
end

% Apply manual CFAR
num_train = 12;
num_guard = 2;
alpha = 10;
detection = cfar_detect(mags, num_train, num_guard, alpha); 



% Plot
x_axis = t.* c;
plot(x_axis, mags, 'b');  % Impulse response
hold on;

% Mark detections
plot(x_axis(detection == 1) , mags(detection == 1), 'ro', 'MarkerFaceColor', 'r');
hold off;

xlabel('Distance from target (m)');
ylabel('Amplitude');
title(sprintf('Impulse Response of multiple targets'));
    




