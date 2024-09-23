function [x,y,u] = ipropFR(U,lambda,z,squaresize)
% Returns Fresnel diffraction pattern of initial field u1 in distance z.
    
    [M, ~] = size(U); % number of samples along single dimension
    size_vector = linspace(-0.5, 0.5, M) * (M-1);
    samp_freq = 1 / squaresize; % sampling frequency
    fx = size_vector / (M-1) * samp_freq;  % frequency vector

    [x, y] = meshgrid(size_vector * squaresize); % spatial coordinates in input plane
    [x_freq, y_freq] = meshgrid(fx * lambda * z); % frequency coordinates in the output plane

    H = exp(-1i * pi / (lambda * z) * (x.^2 + y.^2));  % Inverse quadratic phase factor
    c = exp(1i * 2 * pi / lambda * z) * exp(-1i * pi * (x_freq.^2 + y_freq.^2) / (lambda * z)) / (1i * lambda * z);
    
    u = (ifft2(ifftshift(U./c))) ./H; % Recovered input field

end 
