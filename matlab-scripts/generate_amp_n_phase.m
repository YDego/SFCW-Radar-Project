function [amp, phase] = generate_amp_n_phase(f, sigma, mu, epsilon, r)
    
    amp = zeros(size(f));
    phase = zeros(size(f));

    for i = 1:numel(f)
        temp = 1 + sqrt(1 + (sigma / (2*pi*f(i) * epsilon)) ^ 2);

        alpha = sigma * sqrt((mu / (2 * epsilon)) / temp);
        beta = 2*pi*f(i) * sqrt((mu * epsilon / 2) * temp);
        
        amp(i) = exp(-alpha * 2 * r);
        phase(i) = -beta * 2 * r;
    end

end
