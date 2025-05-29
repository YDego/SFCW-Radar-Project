%% test connection
disp(findPlutoRadio)
% configurePlutoRadio() 


%%

fs = 100e3; % 100 kHz

txPluto = sdrtx('Pluto','RadioID','usb:0','CenterFrequency',1e9, ...
               'BasebandSampleRate',fs,'ChannelMapping',1);
rxPluto = sdrrx('Pluto','RadioID','usb:0','CenterFrequency',1e9, ...
               'BasebandSampleRate',fs,'SamplesPerFrame', 4096,'ChannelMapping',1);
modObj = comm.DPSKModulator('BitInput',true);


%%

% data = exp(1j*2*pi*10e3*(0:1023)/fs).';
t = (0:1023)/fs;
f = 30e3;
sine_1 = (sin(2*pi*f*t)).';         % 10 kHz tone
txData = complex(sine_1);

transmitRepeat(txPluto, txData);


%%

numFrames = 100;
phases = zeros(numFrames,1);
f = 10e3; % tone frequency
N = rxPluto.SamplesPerFrame;
fs = rxPluto.BasebandSampleRate;

% Frequency axis for locating tone bin
freqs = fs * (0:N-1)/N;
[~, binIndex] = min(abs(freqs - f));  % find closest bin to 10kHz


% Collect and analyze frames
for i = 1:numFrames
    rxData = rxPluto();
    spectrum = fft(rxData);
    phase = angle(spectrum(binIndex));
    phases(i) = phase;
end

% Plot phase over time
figure;
plot(phases, '-o')
xlabel('Frame Index')
ylabel('Phase [rad]')
title(sprintf('Phase at %.0f Hz over Time', f))
grid on

figure;
histogram(phases)


%% Cleanup
release(txPluto);
release(rxPluto);
