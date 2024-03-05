function [FX,FY,u2] = propagationFR_analytic(u1, lambda, z,squaresize)
    dx = squaresize;
    [M, N] = size(u1); 
    
    % Spatial frequencies
    fx = (-N/2:N/2-1) / (N*dx);
    fy = (-M/2:M/2-1) / (M*dx);
    [FX, FY] = meshgrid(fx, fy);
    
    % Quadratic phase factor 
    H = (exp(-1i * pi * lambda * z * (FX.^2 + FY.^2)));
    u2 = fftshift(fft2((u1 .* H)));
    u2 = u2./max(max(u2)).*max(max(u1));
end