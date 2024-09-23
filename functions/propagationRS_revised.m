function [x,y,u2] = propagationRS_revised(u1, lambda, z, squaresize)
    
    [M,~] = size(u1); % number of samples along single dimension
    size_vector = linspace(-0.5, 0.5, M) * (M-1);
    samp_freq = 1 / squaresize; % sampling frequency
    fx = size_vector / (M-1) * samp_freq;  % frequency vector

    [x, y] = meshgrid(size_vector * squaresize); % spatial coordinates in input plane
    [x_freq, y_freq] = meshgrid(fx); % frequency coordinates in the output plane
    k = 2*pi/lambda;

    % Free space transfer function
    % ifftshift creates correct frequencies positions
    H = ifftshift((exp(-1i * k * z * sqrt(1 - lambda^2 * (x_freq.^2 + y_freq.^2)))));
    H(isnan(H)) = 0; 
    
    U1 = fft2(u1);
    u2 = ifft2(U1 .* H);

end