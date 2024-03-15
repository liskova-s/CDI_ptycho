function [FX,FY,u2] = propagationFR_analytic(u1, lambda, z,squaresize)
    [M,~] = size(u1); 
    dx1 = squaresize;
    L1 = M*dx1; 
   
    L2 = lambda*z/dx1; 
    dx2 = lambda*z/L1;
    
    x2 = -L2/2:dx2:L2/2-dx2; 
    [FX, FY] = meshgrid(x2, x2);

    % Quadratic phase factor 
    H = (exp(-1i * pi * lambda * z * (FX.^2 + FY.^2)));
    u2 = fftshift(fft2((u1 .* H)));
    u2 = u2./max(max(u2)).*max(max(u1));
end
