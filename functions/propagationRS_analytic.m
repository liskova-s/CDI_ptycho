function [FX,FY,u2] = propagationRS_analytic(U, lambda, z, squaresize)
    k = 2*pi/lambda;
    dx = squaresize;
    [M, N] = size(U); 
    
    % Spatial frequencies
    fx = (-N/2:N/2-1) / (N*dx);
    fy = (-M/2:M/2-1) / (M*dx);
    [FX, FY] = meshgrid(fx, fy);
    
    % Quadratic phase factor 
    H = (exp(-1i * k * z * sqrt(1 - lambda^2 * (FX.^2 + FY.^2))));
    H(isnan(H)) = 0; 
    
    u2 = fftshift(fft2((U .* H)));
    u2 = u2./max(max(u2)).*max(max(U));
end