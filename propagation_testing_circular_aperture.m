% This script performs a visual test of propagation functions on a circular
% aperture.

% Tested functions modelling free-space propagation:
% - propagationFF()             Fraunhoffer propagation model
% - propagationFR()             Fresnel propagation model using impulse
%                               response function
% - propagationFR_analytic()    Fresnel propagation model using analytic
%                               quadratic factor
% - propagationRS()             Rayleigh-Sommerfeld propagation model using
%                               impulse response function
% - propagationRS_analytic()    Rayleigh-Sommerfeld propagation model using
%                               analytic transfer function
%..........................................................................

close all;
clear all;
clc;
addpath('functions/')

%..........................................................................
% GENERATE GRID
% (parameters of detector: 2048 x 2048 px, 1 px = 13,4 um -> cca 2,7 cm)
% -> grid 5000 x 5000: (3 x 3 cm)

gridsize = [5000,5000];
squaresize = 6e-6; 

c = generate_coordinates(gridsize,squaresize);                             

%..........................................................................
% SOURCE SETUP (CIRCULAR APERTURE DIFFRACTION)
tic;
lambda = 555e-9;            % VIS green
source = zeros(gridsize);
R = 2e-4;                   % source radius in m
source(c(:,:,1).^2+c(:,:,2).^2 < R^2) = 1;   % circular footprint: R 

r = c(1,:,1);
k = 2*pi/lambda;
A = 1;
disp(['Source setup - elapsed time: ', num2str(toc), ' sec.']);

%imagesc(source)

%..........................................................................
% FREESPACE PROPAGATION
% (to distance z)
z =  1;                                                     % 1 m distance

tic;
propagatedFresnel1 = propagationFR(source,c,lambda,z);
disp(['Propagation Fresnel, impulse response - elapsed time: ', num2str(toc), ' sec.']);
tic;
[FRX,FRY,propagatedFresnel2] = propagationFR_analytic(source,lambda,z,squaresize);
disp(['Propagation Fresnel, quadratic factor - elapsed time: ', num2str(toc), ' sec.']);
tic;
propagatedRS1 = propagationRS(source,c,lambda,z);
disp(['Propagation Rayleigh-Sommerfield, impulse response - elapsed time: ', num2str(toc), ' sec.']);
tic;
[RSX,RSY,propagatedRS2] = propagationRS_analytic(source, lambda, z, squaresize);
disp(['Propagation Rayleigh-Sommerfield, transfer function - elapsed time: ', num2str(toc), ' sec.']);
tic;
[FFX,FFY,propagatedFraun] = propagationFF(source,lambda,z,squaresize);
disp(['Propagation Fraunhoffer - elapsed time: ', num2str(toc), ' sec.']);

%..........................................................................
% VALIDATION

U_theory = exp(1i*k*z).*exp(1i*k.*r.^2./(2.*z))*A./(z.*1i.*lambda).*2.*(besselj(1,k.*R.*r./z)./(k.*R.*r./z)); 
U_theory = abs(U_theory)/max(abs(U_theory)); % normalization

m1 = z*tan(1.22*lambda/(2*R));

%..........................................................................
% VISUALISATION

figure()
sgtitle("Diffraction - circular aperture")
r = c(1,:,1);

    subplot 261
    imagesc(c(1,:,1),c(:,1,2),abs(source))
    title("Source aperture")

% RS propagation pattern¨
    % RS impulse response FFT
    subplot 262
    imagesc(c(1,:,1),c(:,1,2),abs(propagatedRS1));
    title("RS propagation - impulse response");xlabel("[m]");ylabel("[m]");

    % RS transfer function
    propagatedRS_resc = rescale_interpol(propagatedRS2,FFX,gridsize,c);

    subplot 263
    imagesc(c(1,:,1),c(:,1,2),abs(propagatedRS_resc))
    title("RS propagation - transfer function ");xlabel("[m]");ylabel("[m]");

% RS pattern validation 
    % impulse response validation
    subplot 268
    hold on
    scatter(abs(propagatedRS1(:,round(gridsize(1)/2+1))),r,10,DisplayName='Simulated amplitude profile',  MarkerEdgeColor= "#D95319")
    yline(m1, LineStyle="-."); yline(-m1, Linestyle = "-."); xlabel("[m]"); ylabel("Amplitude");
    title("Rayleigh-Sommerfield propagation")

    % transfer function validation
    subplot 269
    hold on
    scatter(abs(propagatedRS_resc(:,round(gridsize(1)/2+1))),r,10,DisplayName='Simulated amplitude profile', MarkerEdgeColor = "#D95319")
    yline(m1, LineStyle="-."); yline(-m1, Linestyle = "-."); xlabel("[m]"); ylabel("Amplitude");
    title("Rayleigh-Sommerfield propagation")
   

% Fresnel pattern propagation
    % Fresnel impulse response
    subplot 264
    imagesc(c(1,:,1),c(:,1,2),abs(propagatedFresnel1))
    title("Fresnel propagation - impulse response");xlabel("[m]");ylabel("[m]");

    % Fresnel quadratic factor
    propagatedFr_resc = rescale_interpol(propagatedFresnel2,FFX,gridsize,c);

    subplot 265
    imagesc(c(1,:,1),c(:,1,2),abs(propagatedFr_resc))
    %imagesc(abs(propagatedFresnel2(start_index:end_index,start_index:end_index)))
    title("Fresnel propagation - quadratic factor");xlabel("[m]");ylabel("[m]");

% Fresnel pattern validation
    % impulse response
    subplot(2,6,10)
    hold on
    scatter(abs(propagatedFresnel1(:,round(gridsize(1)/2+1))),r,10, DisplayName='Simulated amplitude profile',  MarkerEdgeColor = "#D95319")
    yline(m1, LineStyle="-."); yline(-m1, Linestyle = "-."); xlabel("[m]"); ylabel("Amplitude");
    title("Fresnel propagation")

    % quadratic factor
    subplot(2,6,11)
    hold on
    scatter(abs(propagatedFr_resc(:,round(gridsize(1)/2)+1)),r,10, DisplayName='Simulated amplitude profile', MarkerEdgeColor = "#D95319")
    yline(m1, LineStyle="-."); yline(-m1, Linestyle = "-."); xlabel("[m]"); ylabel("Amplitude");
    title("Fresnel propagation");

% Fraunhoffer pattern
    propagatedFraun_resc = rescale_interpol(propagatedFraun,FFX,gridsize,c);

    subplot 266
    imagesc(c(1,:,1),c(:,1,2),abs(propagatedFraun_resc))
    title("Fraunhofer propagation");xlabel("[m]");ylabel("[m]");

% Fraunhoffer pattern validation
    subplot(2,6,12)
    hold on
    scatter(abs(propagatedFraun_resc(:,round(gridsize(1)/2)+1)),r,10, DisplayName='Simulated amplitude profile', MarkerEdgeColor = "#D95319")
    yline(m1, LineStyle="-."); yline(-m1, Linestyle = "-."); xlabel("[m]"); ylabel("Amplitude");
    title("Fraunhofer propagation")


