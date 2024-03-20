% This script performs visualisation of a diffraction pattern dependence on Fresnel number.

% Tested functions modelling free-space propagation:
% - propagationFF()             Fraunhoffer propagation model
% - propagationFR_analytic()    Fresnel propagation model using analytic
%                               quadratic factor
% - propagationRS_analytic()    Rayleigh-S ommerfeld propagation model using
%                               analytic transfer function
%..........................................................................

close all;
clear;
clc;
addpath('../functions/')

fresnel_list = [1e-3,1e-1,1, 1e1,1e3];

%..........................................................................
% GENERATE GRID
% (parameters of detector: 2048 x 2048 px, 1 px = 13,4 um -> cca 2,7 cm)
% -> grid 5000 x 5000: (3 x 3 cm)

gridsize = [5000,5000];
squaresize = 6e-6; 

c = generate_coordinates(gridsize,squaresize);   

propagatedRS = zeros([gridsize(1),gridsize(2),length(fresnel_list)]);
propagatedFresnel = zeros([gridsize(1),gridsize(2),length(fresnel_list)]);
propagatedFraun = zeros([gridsize(1),gridsize(2),length(fresnel_list)]);
FRXY = zeros([gridsize(1),gridsize(2),2,length(fresnel_list)]);
FFXY = zeros([gridsize(1),gridsize(2),2, length(fresnel_list)]);
RSXY = zeros([gridsize(1),gridsize(2),2, length(fresnel_list)]);

% SOURCE SETUP (CIRCULAR APERTURE DIFFRACTION)
    tic;
    R = 2e-4;           % source radius in m
    lambda = 555e-9;            % VIS green 
    source = zeros(gridsize);
    source(c(:,:,1).^2+c(:,:,2).^2 < R^2) = 1;   % circular footprint: R 
    
    k = 2*pi/lambda;
    A = 1;
    disp(['Source setup - elapsed time: ', num2str(toc), ' sec.']);


% Loop through tested Fresnel numbers:
for i=1:length(fresnel_list)

    fprintf("\nEntering loop for i = %d/%d\n",i,length(fresnel_list));
    F = fresnel_list(i);

    % FREESPACE PROPAGATION
    % F = a^2/(lambda*z) -> z = a^2/(F*lambda)
    z =  R^2/(F*lambda);    


    tic;
    [RSXY(:,:,1,i),RSXY(:,:,2,i),propagatedRS(:,:,i)] = propagationRS_analytic(source,lambda,z,squaresize);
    disp(['Propagation Fresnel, quadratic factor - elapsed time: ', num2str(toc), ' sec.']);
    tic;
    [FRXY(:,:,1,i),FRXY(:,:,2,i),propagatedFresnel(:,:,i)] = propagationFR_analytic(source,lambda,z,squaresize);
    disp(['Propagation Fresnel, quadratic factor - elapsed time: ', num2str(toc), ' sec.']);
    tic;
    [FFXY(:,:,1,i),FFXY(:,:,2,i),propagatedFraun(:,:,i)] = propagationFF(source,lambda,z,squaresize);
    disp(['Propagation Fraunhoffer - elapsed time: ', num2str(toc), ' sec.']);

end

%..........................................................................
% CALCULATE ERROR

disp("\nSubtractions in  progress ..")
diff_Fresnel = zeros([size(propagatedFraun),length(fresnel_list)]);
diff_Fraunhofer = zeros([size(propagatedFraun),length(fresnel_list)]);
mse_Fresnel = zeros(length(fresnel_list));
mse_Fraunhofer = zeros(length(fresnel_list));
diff_Fresnel_amp = zeros([size(propagatedFraun),length(fresnel_list)]);
diff_Fraunhofer_amp = zeros([size(propagatedFraun),length(fresnel_list)]);
mse_Fresnel_amp = zeros(length(fresnel_list));
mse_Fraunhofer_amp = zeros(length(fresnel_list));

for i = 1:length(fresnel_list)
    fprintf("%d/%d\n",i,length(fresnel_list))
    fprintf("  complex field MSE\n")
    diff_Fresnel(:,:,i) = propagatedFresnel(:,:,i)-propagatedRS(:,:,i);
    diff_Fraunhofer(:,:,i) = propagatedFraun(:,:,i)-propagatedRS(:,:,i);

    mse_Fresnel(i) = sum(sum(abs(diff_Fresnel(:,:,i).^2)))/(length(diff_Fresnel)^2);
    mse_Fraunhofer(i) = sum(sum(abs(diff_Fraunhofer(:,:,i).^2)))/(length(diff_Fraunhofer)^2);

    fprintf("  amplitude only MSE\n")
    diff_Fresnel_amp(:,:,i) = abs(propagatedFresnel(:,:,i))-abs(propagatedRS(:,:,i));
    diff_Fraunhofer_amp(:,:,i) = abs(propagatedFraun(:,:,i))-abs(propagatedRS(:,:,i));

    mse_Fresnel_amp(i) = sum(sum(abs(diff_Fresnel_amp(:,:,i).^2)))/(length(diff_Fresnel)^2);
    mse_Fraunhofer_amp(i) = sum(sum(abs(diff_Fraunhofer_amp(:,:,i).^2)))/(length(diff_Fraunhofer)^2);


end

%global_max = max(max(abs(diff_Fresnel),[],"all"),max(abs(diff_Fraunhofer),[],"all"));
%global_min = min(min(abs(diff_Fresnel),[],"all"),min(abs(diff_Fraunhofer),[],"all"));

%..........................................................................
% VISUALISATION
sdisp = @(x) strtrim(evalc(sprintf('disp(%g)', x))); % helper display function


disp("Preparing visualisation Fresnel pattern...");

figure()
sgtitle("Complex pattern subtraction: Fresnel - RS")
for i = 1:length(fresnel_list)
    fprintf("___ 1/2: %d/%d\n",i,length(fresnel_list));
    ax1 = subplot(4,length(fresnel_list),i);
    imagesc(abs(propagatedFresnel(:,:,i)))
    title(sprintf("Fresnel pattern:\nFresnel number: %d\nz = %s m",fresnel_list(i),sdisp(R^2/(fresnel_list(i)*lambda))));
    colormap(ax1,"parula");colorbar;

    ax2 = subplot(4,length(fresnel_list),i+length(fresnel_list));
    imagesc(abs(propagatedRS(:,:,i)))
    colormap(ax2,"parula"); colorbar;
    title("Rayleigh-Sommerfeld pattern")

    ax3 = subplot(4,length(fresnel_list),i+2*length(fresnel_list));
    imagesc(abs(diff_Fresnel(:,:,i)))
    colormap(ax3,"lines");colorbar;
    title(sprintf("Difference complex (FR-RS)\nMSE: %s",sdisp(mse_Fresnel(i))));
  
    ax4 = subplot(4,length(fresnel_list),i+3*length(fresnel_list));
    imagesc(abs(diff_Fresnel(:,:,i)))
    colormap(ax4,"turbo"); colorbar;
end

figure()
sgtitle("Amplitude only subtraction: Fresnel - RS")
for i = 1:length(fresnel_list)
    fprintf("___ 2/2: %d/%d\n",i,length(fresnel_list));
    ax1 = subplot(4,length(fresnel_list),i);
    imagesc(abs(propagatedFresnel(:,:,i)))
    title(sprintf("Fresnel pattern:\nFresnel number: %d\nz = %s m",fresnel_list(i),sdisp(R^2/(fresnel_list(i)*lambda))));
    colormap(ax1,"parula");colorbar;

    ax2 = subplot(4,length(fresnel_list),i+length(fresnel_list));
    imagesc(abs(propagatedRS(:,:,i)))
    colormap(ax2,"parula"); colorbar;
    title("Rayleigh-Sommerfeld pattern")

    ax3 = subplot(4,length(fresnel_list),i+2*length(fresnel_list));
    imagesc(abs(diff_Fresnel_amp(:,:,i)))
    colormap(ax3,"lines");colorbar;
    title(sprintf("Difference amplitude (FR-RS)\nMSE: %s",sdisp(mse_Fresnel_amp(i))));

    ax4 = subplot(4,length(fresnel_list),i+3*length(fresnel_list));
    imagesc(abs(diff_Fresnel_amp(:,:,i)))
    colormap(ax4,"turbo"); colorbar;
end
%__________________________________________________________________________
disp("Preparing visualisation Fraunhofer...");

figure()
sgtitle("Complex pattern subtraction: Fraunhofer - RS")
for i = 1:length(fresnel_list)
    fprintf("___ 1/2: %d/%d\n",i,length(fresnel_list));
    ax1 = subplot(4,length(fresnel_list),i);
    imagesc(abs(propagatedFraun(:,:,i)))
    title(sprintf("Fraunhofer pattern:\nFresnel number: %d\nz = %s m",fresnel_list(i),sdisp(R^2/(fresnel_list(i)*lambda))));
    colormap(ax1,"parula");colorbar;

    ax2 = subplot(4,length(fresnel_list),i+length(fresnel_list));
    imagesc(abs(propagatedRS(:,:,i)))
    colormap(ax2,"parula"); colorbar;
    title("Rayleigh-Sommerfeld pattern")

    ax3 = subplot(4,length(fresnel_list),i+2*length(fresnel_list));
    imagesc(abs(diff_Fraunhofer(:,:,i)))
    colormap(ax3,"lines");colorbar;
    title(sprintf("Difference complex (FR-RS)\nMSE: %s",sdisp(mse_Fraunhofer(i))));

    ax4 = subplot(4,length(fresnel_list),i+3*length(fresnel_list));
    imagesc(abs(diff_Fraunhofer(:,:,i)))
    colormap(ax4,"turbo"); colorbar;
end

figure()
sgtitle("Amplitude only subtraction: Fraunhofer - RS")
for i = 1:length(fresnel_list)
    fprintf("___ 2/2: %d/%d\n",i,length(fresnel_list));
    ax1 = subplot(4,length(fresnel_list),i);
    imagesc(abs(propagatedFraun(:,:,i)))
    colormap(ax1,"parula");colorbar;
    title(sprintf("Fraunhofer pattern:\nFresnel number: %d\nz = %s m",fresnel_list(i),sdisp(R^2/(fresnel_list(i)*lambda))));
    
    ax2 = subplot(4,length(fresnel_list),i+length(fresnel_list));
    imagesc(abs(propagatedRS(:,:,i)))
    colormap(ax2,"parula"); colorbar;
    title("Rayleigh-Sommerfeld pattern")

    ax3 = subplot(4,length(fresnel_list),i+2*length(fresnel_list));
    imagesc(abs(diff_Fraunhofer_amp(:,:,i)))
    colormap(ax3,"lines");colorbar;
    title(sprintf("Difference amplitude (FR-RS)\nMSE: %s",sdisp(mse_Fraunhofer_amp(i))));

    ax4 = subplot(4,length(fresnel_list),i+3*length(fresnel_list));
    imagesc(abs(diff_Fraunhofer_amp(:,:,i)))
    colormap(ax4,"turbo"); colorbar;
end
%_______________________________
