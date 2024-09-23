% Script that calculates and displays diffraction profile for a set of
% Fresnel numbers and a pair of grid parameters. Including MSE calculation (RS pattern as reference). 
% This script enables fine-tuning of the combination of Fresnel number,
% gridsize and squaresize parameters to eliminate numerical ocillations
% inside the pattern.

close all;
clear;
clc;
addpath('../functions/')

% FIXED PARAMETERS
lambda = 555e-9;  % m 
R = 500e-6; % source radius in m
k = 2*pi/lambda;

% Fresnel number
fresnel_list = [0.945,0.95]; % <- 
z_list = R^2./(fresnel_list.*lambda); 

% GRID parameters
gridsize = [8001,8001]; % <- 
squaresize = 0.00008; % <- 

%_________________________________________________
c = generate_coordinates(gridsize,squaresize);   
source = zeros(gridsize);
source(c(:,:,1).^2+c(:,:,2).^2 < R^2) = 1;  

tic;
mse = propagated_mse(gridsize,fresnel_list,z_list,source,lambda,squaresize,length(fresnel_list));
toc

% FUNCTIONS________________________________________________________________
function out = scale_energy(propagated_intensity, xF)
    energy = trapz(xF, propagated_intensity); %trapezoidal integration
    out = propagated_intensity / energy;
end

function mse = getMSE(vector, reference)
    mse = sum((abs(vector) - abs(reference)).^2)./sum(abs(reference).^2);
end

function mse_vector = propagated_mse(gridsize,fresnel_list,z_list,source,lambda,squaresize,N)
    fprintf("\nSTARTING PROPAGATION\n");
    mse_vector = zeros([N,3]); %(:,1) for Fresnel, (:,2) for Fraunhofer
    % for each fresnel number perform whole MSE calculation
    
    for i=1:N
        % PROPAGATION
        fprintf("Loop %d/%d\n",i,N)
        z = z_list(i);
        [x,~,pRS] = propagationRS_revised(source,lambda,z,squaresize);
        xRS= x(1,:);
        fprintf("             RS propagation completed.\n")
        [x,~,pFresnel] = propagationFR_revised(source,lambda,z,squaresize);
        xFR = x(1,:);
        fprintf("             FRESNEL propagation completed.\n")
        [x,~,pFraun] = propagationFF(source,lambda,z,squaresize);
        xFF = x(1,:);
        fprintf("             Fraunhofer propagation completed.\n")
        
        % CUT 1D PROFILES
        cut = ceil(gridsize(1)/2);
        prop_RS=squeeze(abs(pRS(cut,:).^2));
        prop_Fresnel=squeeze(abs(pFresnel(cut,:).^2));
        prop_Fraunhofer=squeeze(abs(pFraun(cut,:).^2));

        % NORMALIZATION
        profiles_RS = scale_energy(prop_RS,xRS);
        profiles_Fresnel = scale_energy(prop_Fresnel,xFR);
        profiles_Fraunhofer = scale_energy(prop_Fraunhofer,xFF);
        
        % VISUALISATION
        figure()
        plot(xRS,profiles_RS,"DisplayName","RS")
        hold on
        plot(xFR,profiles_Fresnel,"DisplayName","FR")
        plot(xFF,profiles_Fraunhofer,"DisplayName","FF")
        title(sprintf("Fresnel number: %f",fresnel_list(i)))
        legend()

        % MSE CALCULATION:
    
        % interpolation 
        if xFF(1) > xRS(1)
            % xRS coarser
            intensityFF = interp1(xFF, profiles_Fraunhofer, xRS, 'linear', 'extrap'); 
            intensityFR = interp1(xFR, profiles_Fresnel, xRS, 'linear', 'extrap');
            intensityRS = profiles_RS;
        else
            intensityFF = profiles_Fraunhofer;
            intensityFR = profiles_Fresnel;
            intensityRS = interp1(xRS, profiles_RS, xFR, 'linear', 'extrap');
        end
       
        % defining the area of interestfor MSE calculation:
        % given by spatial extent of FR and FF patterns (as it remains equal and is given by Fresnel number)
        startIdx = find(xRS >= min(xFF), 1, 'first');
        endIdx = find(xRS <= max(xFF), 1, 'last');
    
        % Crop vectors
        croppedIntensityFF_RS = intensityFF(startIdx:endIdx);
        croppedIntensityFR_RS = intensityFR(startIdx:endIdx);
        croppedIntensityRS = intensityRS(startIdx:endIdx);

        mse_vector(i,1) = fresnel_list(i);
        mse_vector(i,2) = getMSE(croppedIntensityFR_RS,croppedIntensityRS);
        mse_vector(i,3) = getMSE(croppedIntensityFF_RS,croppedIntensityRS);
        fprintf("             MSE %.10f %.10f\n",mse_vector(i,2),mse_vector(i,3))
        
        %save("mse_2-5_new.mat","mse_vector")
    end
end



