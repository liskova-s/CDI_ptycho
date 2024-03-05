function [X2,Y2,u2] = propagationFF(u1, lambda, z,squaresize)
    [M,~] = size(u1); 
    dx1 = squaresize;
    L1 = M*dx1; 
   
    L2 = lambda*z/dx1; 
    dx2 = lambda*z/L1;
    
    x2 = -L2/2:dx2:L2/2-dx2; 
    [X2, Y2] = meshgrid(x2, x2);
    
    c = exp(1i*2*pi/lambda*z) * exp(1i*pi*(X2.^2 + Y2.^2)/(lambda*z)) / (1i*lambda*z);
    u2 = c .* fftshift(fft2(fftshift(u1))) * dx1^2;
    
    u2 = u2 / max(max(abs(u2)));
end