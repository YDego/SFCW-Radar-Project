function [amp, phase] = generate_amp_n_phase(w, sigma, mu, epsilon, r)
    
    sub_function = generate_sub_function(w, epsilon, sigma);

    alpha = sigma * sqrt((mu / 2 / epsilon) ./ sub_function);
    beta = w .* sqrt((mu * epsilon / 2) * sub_function);

    phase = -beta .* (2 * r);
    amp = exp(-alpha .*  2 * r);

end
