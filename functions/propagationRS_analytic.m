function [FX,FY,u2] = propagationRS_analytic(U, lambda, z, squaresize)
    [M,~] = size(U); 
    dx1 = squaresize;
    L1 = M*dx1; 
   
    L2 = lambda*z/dx1; 
    dx2 = lambda*z/L1;
    
    x2 = -L2/2:dx2:L2/2-dx2; 
    [FX, FY] = meshgrid(x2, x2);
    k = 2*pi/lambda;
    % Quadratic phase factor 
    H = (exp(-1i * k * z * sqrt(1 - lambda^2 * (FX.^2 + FY.^2))));
    H(isnan(H)) = 0; 
    
    u2 = fftshift(fft2((U .* H)));
    u2 = u2./max(max(u2)).*max(max(U));
    u2 = u2 / max(max(abs(u2)));
end
