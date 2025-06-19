function detection = cfar_detect(signal, num_train, num_guard, alpha)
    %CFAR_DETECT Apply 1D CA-CFAR on the input signal
    %   detection = cfar_detect(signal, num_train, num_guard, alpha)
    %
    %   Inputs:
    %       signal     - input signal (e.g., abs(impulse response))
    %       num_train  - number of training cells on each side
    %       num_guard  - number of guard cells on each side
    %       alpha      - threshold scaling factor (higher → fewer false alarms)
    %
    %   Output:
    %       detection  - binary vector: 1 if target detected, 0 otherwise
    
    N = length(signal);
    detection = zeros(1, N);
    
    for i = num_train + num_guard + 1 : N - num_train - num_guard
        % Define training cells (excluding guard cells)
        training_left = signal(i - num_guard - num_train : i - num_guard - 1);
        training_right = signal(i + num_guard + 1 : i + num_guard + num_train);
        
        % Estimate noise level from training cells
        noise_level = mean([training_left, training_right]);
        
        % CFAR threshold
        threshold = alpha * noise_level;
        
        % If signal exceeds threshold, declare detection
        if signal(i) > threshold
            detection(i) = 1;
        end
    end

end
