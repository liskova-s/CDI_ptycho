function [x_freq,y_freq,u2] = propFR(u1,lambda,z,squaresize)
% Returns Fresnel diffraction pattern of initial field u1 in distance z.
    
    [M, ~] = size(u1); % number of samples along single dimension
    size_vector = linspace(-0.5, 0.5, M) * (M-1);
    samp_freq = 1 / squaresize; % sampling frequency
    fx = size_vector / (M-1) * samp_freq;  % frequency vector

    [x, y] = meshgrid(size_vector * squaresize); % spatial coordinates in input plane
    [x_freq, y_freq] = meshgrid(fx * lambda * z); % frequency coordinates in the output plane

    H = exp(-1i * pi / (lambda * z) * (x.^2 + y.^2)); % Quadratic phase factor 
    c = exp(1i * 2 * pi / lambda * z) * exp(-1i * pi * (x_freq.^2 + y_freq.^2) / (lambda * z)) / (1i * lambda * z);
    u2 = c .* fftshift(fft2((u1 .* H))); % Output field

end 

