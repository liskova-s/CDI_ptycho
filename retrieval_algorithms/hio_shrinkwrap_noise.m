% Simulating and measuring noise effect on hio + shrinkwrap reconstruction.

%close all;
clear;
clc;

addpath('../../functions/')

% Far-field zone: N_f<<1
lambda = 555e-9; 
squaresize = 1e-6;
z = 10;

% get object and diffraction pattern
imageorig = abs(cell2mat(struct2cell(load("logo.mat")))); % nophase fjfi logo
imageorig  = imbinarize(imresize(imageorig , 0.2));
padding = 723;
image = padarray(imageorig , [padding, padding], 'both');
oversampling_ratio = size(image)/size(imageorig);

[~,~,FT]= propFR(image,lambda,z,squaresize); 
intensity = abs(FT).^2;

% initial random phase
random_phase = 2*pi*rand(size(FT));

% form object domain constraint mask
object_mask = zeros(size(FT));
object_mask(padding:end-padding,padding:end-padding)=1;

snr_list = [300,500]; % specify simulated noise as SNR (multiple simulations at once)
N_iter = 1500;

error_object_dom = zeros([length(snr_list),N_iter]);
error_object_MSE = zeros([length(snr_list),N_iter]);
error_object_SSIM = zeros([length(snr_list),N_iter]);
error_object_NCC = zeros([length(snr_list),N_iter]);
error_fourier_MSE = zeros([length(snr_list),N_iter]);
    
for j = 1:length(snr_list)
    % introduce noise to intensity pattern
    SNR = snr_list(j);
    signal_power = mean(intensity(:).^2);
    noise_power = signal_power / SNR; % noise power is it's std
    
    rng(1, "twister");
    noise = sqrt(noise_power) * randn(size(intensity)); % zero mean, adjusted variance
    noisy_intensity = intensity + noise;
    
    measured_amp = abs(sqrt(noisy_intensity)); 
    measured_object = abs(double(image)); 
   
    % create initial input
    fprintf("Creating initial input.\n")
    G_0 = measured_amp.*exp(1i.*random_phase);    %initial detector field
    [~,~,g_0]= ipropFR(G_0,lambda,z,squaresize); 
    input_g = g_0;
    
    figure()
    % HIO iterations___________________________________________________________
    beta = 0.8;
    shrinkwrap_sigma = 5;
    shrinkwrap_freq = 20; % update support every {shrinkwrap_freq}-th iteration
    end_sigma = 2.5;
    delta_sigma = (shrinkwrap_sigma - end_sigma)/(N_iter/shrinkwrap_freq - 1);
    
    count = 1;
    figure()
    shrinkwrap_count = 1;
    fprintf("Entering iterative loop.\n")
  
    for k = 1:N_iter
        fprintf("iteration %d/%d\n",k,N_iter)
        %1) FT of input
        [~,~,G_k]= propFR(input_g,lambda,z,squaresize); 
        error_fourier_MSE(k) = error_intensity(abs(G_k).^2, abs(measured_amp).^2);
        
        %2) replace amplitude with original one
        G_kk = measured_amp.*exp(1i.*angle(G_k));
    
        % 3) IFFT back to object domain
        [~,~,output_g]= ipropFR(G_kk,lambda,z,squaresize); 
        error_object_dom(j,k) = error_object_domain(output_g,object_mask);
        error_object_MSE(j,k) = error_MSE(measured_object, abs(output_g));
        error_object_SSIM(j,k)= error_SSIM(measured_object, abs(output_g));
        error_object_NCC(j,k) = error_NCC(measured_object, abs(output_g));
        
        % SHRINKWRAP support update
        if mod(k,shrinkwrap_freq)==0
            % calculate gaussian
            N = size(object_mask,1);
            [X, Y] = meshgrid(-N/2+1:N/2, -N/2+1:N/2);
            sigma = shrinkwrap_sigma - shrinkwrap_count*delta_sigma;
            G = 1 * exp(-(X.^2 + Y.^2) / (2 * sigma^2));
    
            %[~,~,ffG]= propFR(G,lambda,z,squaresize); 
            %[~,~,iffg]= propFR(input_g,lambda,z,squaresize); 
            %[~,~,ob_mask]= ipropFR(ffG.*iffg,lambda,z,squaresize); 
            ob_mask = abs(ifftshift(ifft2(fft2(fftshift(G)).*fft2(fftshift(abs(input_g))))));
            ob_mask = rescale(abs(ob_mask), 0, 1);
            
            temp_mask = zeros(size(ob_mask));
            temp_mask(ob_mask>0.12) = 1;
            object_mask = object_mask.*temp_mask;
            shrinkwrap_count = shrinkwrap_count + 1;
        end
    
    % constraints in the object plane
        % 4) positivity + object constraints: form new input
        % object constaint
        new_g = output_g;
        mask_violation = (abs(output_g) <= 0) | (object_mask == 0);
        new_g(mask_violation) = input_g(mask_violation) - beta * output_g(mask_violation);
        new_g(~mask_violation) = output_g(~mask_violation);
        input_g = new_g;
        
        % visualisation
        if mod(k,round(N_iter/5))==0
            subplot(1,5,count)
            count=count+1;
            imagesc(abs(output_g))
            title(sprintf("Iteration number %d",k))
        end
    end
sgtitle("HIO reconstruction")
figure()
end

figure()
plot(error_object_dom(1,:),"DisplayName","Fienup's error")
hold on
plot(1 - error_object_SSIM(1,:),"DisplayName","1 - SSIM")
plot(1 - error_object_NCC(1,:),"DisplayName","1 - NCC")
title("Error evolution in object domain")
xlabel("Iterations")
ylabel("Error")
grid on
legend()

figure()
plot(error_object_MSE(1,:),"DisplayName","MSE")
xlabel("Iterations")
ylabel("MSE")
grid on
title("MSE evolution")

figure()
plot(error_fourier_MSE(1,:),"DisplayName","MSE fourier dom")
xlabel("Iterations")
ylabel("MSE")
grid on
title("MSE evolution fourier domain")

function error = error_object_domain(g_k,mask)
    % RMS
    % G...calculated FT, F...measured pattern

    num = abs(g_k.*~mask).^2;
    denom = abs(g_k).^2;
    error = sqrt(sum(num, 'all') / sum(denom, 'all'));
end

function mse = error_MSE(original, reconstructed)
% compute MSE between current reconstruction and target
    mse = mean((original - reconstructed).^2, 'all');
end

function ssim_val = error_SSIM(original, reconstructed)
% evaluate structural similarity index
    [ssim_val, ~] = ssim(reconstructed, original);
end

function ncc = error_NCC(original, reconstructed)
% evaluate normalized cross-correlation
    %ncc = sum((original(:) .* reconstructed(:))) / (norm(original(:)) * norm(reconstructed(:)));
    % Calculate the means of the original and reconstructed signals
    mean_original = mean(original(:));
    mean_reconstructed = mean(reconstructed(:));

    % Calculate the numerator (sum of product of deviations from mean)
    numerator = sum((original(:) - mean_original) .* (reconstructed(:) - mean_reconstructed));

    % Calculate the denominator (product of the square roots of the sums of the squared deviations)
    denominator = sqrt(sum((original(:) - mean_original).^2) * sum((reconstructed(:) - mean_reconstructed).^2));

    % Calculate normalized cross-correlation
    ncc = numerator / denominator;
end

function int_error = error_intensity(reconstructed,measured)
    % MSE of intensity in fourier domain
     R = (abs(reconstructed)./max(max(abs(reconstructed)))).^2;
     M = (abs(measured)./max(max(measured))).^2;
     int_error = mean((R - M).^2, 'all');
end
