clear
clc

%% Parameters

% constants & system parameters
e_0 = 8.8541e-12;       % [F/m]
mu_0 = 1.2566e-6;       % [N/(A^2)]
e_r = 1.00054;
mu_r = 1;
sigma = 1e-9;              % conductivity [S/m]
epsilon = e_r * e_0;
mu = mu_r * mu_0;
c = 1/sqrt(mu * epsilon);

% freq parameters
BW = 5e9;               % 2[GHz]
fs = 2 * BW;            % sample rate

f_start = 1e9;        % [Hz]
f_end = 3e9;          % [Hz]

% inputs
N = 1024;               % window size [samples]

% time & freq vectors
t = 0:1/fs:1/fs*(N-1);
f = linspace(0, BW, N);
w = 2*pi*f;

% Find indices where freq is within the specified range
cut_indices = (f >= f_start) & (f <= f_end);

% Calculate effective bandwidth and number of used frequency points
B = f_end - f_start;   % [Hz]
n = sum(cut_indices);  % number of frequency bins in band

% Frequency step in used region
delta_f = B / n;

% Range resolution and max range
range_resolution = c / (2 * B);          % [m]
r_max = c / (2 * delta_f);               % [m]
r_min = range_resolution;
%% amplitude & phase
% Calculate amplitude and phase
[amp, phase] = generate_amp_n_phase(w, sigma, mu, epsilon, r_max);

% Plot the original and filtered amplitude vs frequency
figure(1);
sgtitle(sprintf('Wave propagation (vacum), z = %f(m)', 2*r_max))
subplot(2, 1, 1);
plot(f(cut_indices)*1e-9, 20*log(amp(cut_indices)));
xlabel('Frequency (GHz)');
ylabel('Amplitude (dB)');
title('Amplitude vs Frequency');

subplot(2, 1, 2+);
plot(f(cut_indices)*1e-9, phase(cut_indices));
xlabel('Frequency (GHz)');
ylabel('Phase (radians)');
title('Phase vs Frequency');


%% Create a single figure with 3 subplots


% inputs
r_values = linspace(r_min, r_max, 4);   % target distances [m]
r_values = r_values(2:end-1);
figure(2);

for i = 1:length(r_values)
    r = r_values(i);
    
    % Calculate amplitude and phase
    [amp, phase] = generate_amp_n_phase(w, sigma, mu, epsilon, r);
    
    % Compute frequency response and impulse response
    freq_response = amp .* exp(1j * phase);
    freq_response(~cut_indices) = 0;
    imp_response = ifft(freq_response, 'symmetric');
    mag = abs(imp_response);

    % % Plot impulse response in each subplot
    % subplot(length(r_values), 1, i);
    % plot(t*1e6, mag);
    % xlabel('Time (\mu s)');
    % ylabel('Amplitude');
    % title(sprintf('Impulse Response (r = %d m)', r));

    % % Apply manual CFAR
    %     num_train = 20;
    %     num_guard = 4;
    %     alpha = 5;
    %     detection = cfar_detect(mag, num_train, num_guard, alpha); 
    
    % Find index of the center
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

    xlabel('Distance (m)');
    ylabel('Amplitude');
    title(sprintf('Impulse Response (r = %.3f m)', r));
    
end


