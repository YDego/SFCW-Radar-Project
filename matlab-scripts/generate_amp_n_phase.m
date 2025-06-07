function [amp, phase] = generate_amp_n_phase(w, sigma, mu, epsilon, r)
    
    amp = zeros(size(w));
    phase = zeros(size(w));

    for i = 1:numel(w)
        temp = 1 + sqrt(1 + (sigma / (w(i) * epsilon)) ^ 2);

        alpha = sigma * sqrt((mu / (2 * epsilon)) / temp);
        beta = w(i) * sqrt((mu * epsilon / 2) * temp);
        
        amp(i) = exp(-alpha * 2 * r);
        phase(i) = -beta * (2 * r);
    end

end
