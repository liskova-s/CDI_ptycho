function u2 = propagationRS(u1,coordinates, lambda, z)
    k = 2*pi/lambda;
    rsq = (coordinates(:,:,1).^2+coordinates(:,:,2).^2+z.^2);
    % impulse response
    h = z.*exp(1i*k*sqrt(rsq))./(1i.*lambda.*rsq);

    % padding
    mm = size(u1,1) + size(h,1) - 1;
    nn = size(u1,2) + size(h,2) - 1;
    
    % angular spectrum 
    U = fft2(u1, mm, nn); % Fourier transform of the field
    H = fft2(h, mm, nn);  % FFT of impulse response function 
    C = ifft2(U .* H);

    % padding removal
    pad_m = ceil((size(h,1)-1) / 2);
    pad_n = ceil((size(h,2)-1) / 2);
    u2 = C(pad_m+1:size(u1,1)+pad_m, pad_n+1:size(u1,2)+pad_n);
    u2 = u2./max(max(u2)).*max(max(u1));
end