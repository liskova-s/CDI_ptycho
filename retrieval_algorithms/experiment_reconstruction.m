% This script simulates diffraction pattern corresponding to
% parameters of the experiment, considering dynamic range of camera and
% performs hybrid input-output phase retrieval using Fresnel propagation (propFR, ipropFR).

%close all;
clear;
clc;
addpath('../../functions/')

% parameters
lambda = 75/21*1e-9; % [m] 21st harmonic
z = 100e-3; % between 100 - 150 mm
R = 5e-6/2; % object size: 5 - 10 um
% N_f ~ 0.0175 - 0.186

% camera chip: 2024x2024 px
% pixel size 15.6 um
% size of measured intensity pattern: 27.5x27.5 mm

squaresize = 6.46e-9; % to obtain diff pattern of size of chip
% source object has to have 2R/squaresize pixels -> 774x774 to represent 5
% um in real space

imageorig = abs(cell2mat(struct2cell(load("logo.mat")))); % nophase fjfi logo
imageorig  = imbinarize(imresize(imageorig , 0.2)); % scaling to 774x774
% ensure oversampling (sigma ~ 12)
padding = 972; % 2048 x 2048
image = padarray(imageorig , [padding, padding], 'both');
oversampling_ratio = size(image)/size(imageorig)

% fine diffraction pattern - 10000 x 10000
[xFT,~,FT]= propFR(image,lambda,z,squaresize); 

% Forming camera chip output (measured_intensity)
% sample to 2024x2024 (camera chip output) + adjust intensity scale
% 16 bit camera - 65536 gray levels
intensity = abs(FT).^2;
scaled_vector = rescale(intensity, 0, 65000); %16bit (65535) ale hodnota maxima je kolem 45000
intensity  = double(uint16(scaled_vector));

% Form object mask - square of object size
% current matrix size: 2024x2024
object_mask = zeros(size(intensity));
object_mask(padding:end-padding,padding:end-padding)=1;

% add noise
sigma = 20;
noise_mean = 100;
signal_power = mean(intensity(:).^2);
rng(1, "twister");
noise = noise_mean + sigma.*randn(size(intensity)); 
noise_power = (sigma.^2);
SNR = signal_power/noise_power
noisy_intensity = double(intensity) + noise;
scaled_vector = rescale(noisy_intensity, 0, 65000); %16bit (65535) ale hodnota maxima je kolem 45000
noisy_intensity  = double(uint16(scaled_vector));
%noisy_intensity = intensity;

measured_amp = abs(sqrt(double(noisy_intensity)));
measured_object = abs(double(image));

% HIO ______________________________________________________________________
beta = 0.8;
N_iter = 2400;
shrinkwrap_sigma = 7; %5
shrinkwrap_freq = 20; % update support every {shrinkwrap_freq}-th iteration
end_sigma = 2.5;
delta_sigma = (shrinkwrap_sigma - end_sigma)/(N_iter/shrinkwrap_freq - 1);

count = 1;
figure()
shrinkwrap_count = 1;
shrinkwrap_threshold = 0.18; % 0.12, .18

rng(1,"twister");
random_phase = 2*pi*rand(size(measured_amp));
G_0 = measured_amp.*exp(1i.*random_phase);    %initial detector field

[~,~,g_0]= ipropFR(G_0,lambda,z,squaresize); 
input_g = g_0;

fprintf("Entering iterative loop.\n")

error_object_dom = zeros([1,N_iter]);
error_object_MSE = zeros([1,N_iter]);
error_object_SSIM = zeros([1,N_iter]);
error_object_NCC = zeros([1,N_iter]);

for k = 1:N_iter
    fprintf("iteration %d/%d\n",k,N_iter)
    %1) FT of input
    [~,~,G_k]= propFR(input_g,lambda,z,squaresize); 
    
    %2) replace amplitude with original one
    G_kk = measured_amp.*exp(1i.*angle(G_k));

    % 3) IFFT back to object domain
    [~,~,output_g]= ipropFR(G_kk,lambda,z,squaresize); 

    error_object_dom(k) = error_object_domain(output_g,object_mask);
    error_object_MSE(k) = error_MSE(measured_object, abs(output_g));
    error_object_SSIM(k)= error_SSIM(measured_object, abs(output_g));
    error_object_NCC(k) = error_NCC(measured_object, abs(output_g));
    
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
        temp_mask(ob_mask>shrinkwrap_threshold) = 1;
        object_mask = object_mask.*temp_mask;
        shrinkwrap_count = shrinkwrap_count + 1;
    end

    % 4) positivity + object constraints: form new input
    % object constaint
    new_g = output_g;
    mask_violation = (abs(output_g) <= 0) | (object_mask == 0);
    new_g(mask_violation) = input_g(mask_violation) - beta * output_g(mask_violation);
    new_g(~mask_violation) = output_g(~mask_violation);
    
    input_g = new_g;
    
    % visualisation
    if mod(k,round(N_iter/10))==0
        subplot(2,5,count)
        count=count+1;
        imagesc(rot90(rot90(abs( output_g ))))
        title(sprintf("Iteration number %d",k))
    end
end
sgtitle("HIO reconstruction")

figure()
plot(error_object_dom,"DisplayName","Fienup's error")
hold on
plot(1 - error_object_SSIM,"DisplayName","1 - SSIM")
plot(1 - error_object_NCC,"DisplayName","1 - NCC")
title("Error evolution in object domain")
xlabel("Iterations")
ylabel("Error")
grid on
legend()

figure()
plot(error_object_MSE,"DisplayName","MSE")
xlabel("Iterations")
ylabel("MSE")
grid on
title("MSE evolution")

function error = error_object_domain(g_k,mask)
    % RMS
    % G...calculated FT, F...measured pattern

    num = abs(g_k.*~mask).^2;
    denom = abs(g_k).^2;
    error = sqrt(sum(num, 'all') / sum(denom, 'all'));
end

function mse = error_MSE(original, reconstructed)
% compute MSE between current reconstruction and target
    mse1 = mean((original - reconstructed).^2, 'all');
    mse2 = mean((rot90(rot90(original)) - reconstructed).^2, 'all');
    mse = min(mse1,mse2);
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

