function u2 = propagationFR(wave, coordinates, lambda, z)
    k = 2*pi/lambda;
    rsq = (coordinates(:,:,1).^2+coordinates(:,:,2).^2);
    h = exp(1i*k*z)/(1i*lambda*z)*exp(1i*k/(2*z)*rsq);
     % padding
    mm = size(wave,1) + size(h,1) - 1;
    nn = size(wave,2) + size(h,2) - 1;
    
    % convolution via FT
    W = fft2(wave, mm, nn); % Fourier transform of the field
    H = fft2(h, mm, nn);    % FFT of impulse response function 
    C = ifft2(W .* H);

    % padding removal
    pad_m = ceil((size(h,1)-1) / 2);
    pad_n = ceil((size(h,2)-1) / 2);
    result = C(pad_m+1:size(wave,1)+pad_m, pad_n+1:size(wave,2)+pad_n);
    u2 = result./max(max(result)).*max(max(wave));
end
