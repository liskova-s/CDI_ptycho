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
addpath('../functions/')

fresnel_list = [0.1, 0.5, 0.9, 1, 4, 10, 50];

%..........................................................................
% GENERATE GRID
% (parameters of detector: 2048 x 2048 px, 1 px = 13,4 um -> cca 2,7 cm)
% -> grid 5000 x 5000: (3 x 3 cm)

gridsize = [5000,5000];
squaresize = 6e-6; 

c = generate_coordinates(gridsize,squaresize);   
z = 1;

propagatedFresnel1 = zeros([gridsize(1),gridsize(2),length(fresnel_list)]);
propagatedFresnel2 = zeros([gridsize(1),gridsize(2),length(fresnel_list)]);
propagatedFraun = zeros([gridsize(1),gridsize(2),length(fresnel_list)]);
FRXY = zeros([gridsize(1),gridsize(2),2,length(fresnel_list)]);
FFXY = zeros([gridsize(1),gridsize(2),2, length(fresnel_list)]);

% Loop through tested Fresnel numbers:
for i=1:length(fresnel_list)

    fprintf("\nEntering loop for i = %d\n",i);
    F = fresnel_list(i);
    % SOURCE SETUP (CIRCULAR APERTURE DIFFRACTION)
    tic;
    R = 2e-4;           % source radius in m
    lambda = R^2/(F*z); % varying wavelength    
    source = zeros(gridsize);
    source(c(:,:,1).^2+c(:,:,2).^2 < R^2) = 1;   % circular footprint: R 
    
    k = 2*pi/lambda;
    A = 1;
    disp(['Source setup - elapsed time: ', num2str(toc), ' sec.']);

    % FREESPACE PROPAGATION
    % F = a^2/(lambda*z) -> z = a^2/(F*lambda)
    %z =  R^2/(F*lambda);                                             

    tic;
    propagatedFresnel1(:,:,i) = propagationFR(source,c,lambda,z);
    disp(['Propagation Fresnel, impulse response - elapsed time: ', num2str(toc), ' sec.']);
    tic;
    [FRXY(:,:,1,i),FRXY(:,:,2,i),propagatedFresnel2(:,:,i)] = propagationFR_analytic(source,lambda,z,squaresize);
    disp(['Propagation Fresnel, quadratic factor - elapsed time: ', num2str(toc), ' sec.']);
    tic;
    [FFXY(:,:,1,i),FFXY(:,:,2,i),propagatedFraun(:,:,i)] = propagationFF(source,lambda,z,squaresize);
    disp(['Propagation Fraunhoffer - elapsed time: ', num2str(toc), ' sec.']);

end

%..........................................................................
% VISUALISATION
disp("Preparing visualisation 1/2");
figure()
sgtitle("Fraunhofer propagation with changing Fresnel number, circular aperture")
r = c(1,:,1);

for i = 1:length(fresnel_list)
    fprintf("Case %d/%d\n",i,length(fresnel_list));
    subplot(2,length(fresnel_list),i)
    %propagatedFF_resc = rescale_interpol(propagatedFraun(:,:,i),FFXY(:,1,i),gridsize,c);
    propagatedFF_resc= propagatedFraun(:,:,i);
    imagesc(FFXY(1,:,1,i),FFXY(:,1,2,i),abs(propagatedFF_resc))
    title(sprintf("Fresnel number: %.2f\nλ = %.2f nm",fresnel_list(i), 10^9*R^2/(fresnel_list(i)*z)));xlabel("[m]");ylabel("[m]");
    
    subplot(2,length(fresnel_list),i+length(fresnel_list))
    scatter(abs(propagatedFF_resc(:,round(gridsize(1)/2)+1)),FFXY(1,:,1,i),10, MarkerEdgeColor = "#D95319")
    xlabel("[m]"); ylabel("Amplitude");
    title(sprintf("Fresnel number: %.2f\n",fresnel_list(i)));
    
end

disp("Preparing visualisation 2/2");
figure()
sgtitle("Fresnel propagation with changing Fresnel number, circular aperture")
r = c(1,:,1);

for i = 1:length(fresnel_list)
    fprintf("Case %d/%d\n",i,length(fresnel_list));

    subplot(2,length(fresnel_list),i)
    imagesc(FRXY(1,:,1,i),FRXY(:,1,2,i),abs(propagatedFresnel1(:,:,i)))
    xlabel("[m]");ylabel("[m]");
    title(sprintf("Fresnel number: %.2f\nλ = %.3f nm",fresnel_list(i), 10^9*R^2/(fresnel_list(i)*z)));
    
    subplot(2,length(fresnel_list),i+length(fresnel_list))
    scatter(abs(propagatedFresnel1(:,round(gridsize(1)/2)+1,i)),FRXY(1,:,1,i),10, MarkerEdgeColor = "#D95319")
    xlabel("[m]"); ylabel("Amplitude");
    title(sprintf("Fresnel number: %.2f\n",fresnel_list(i)));
    
end
