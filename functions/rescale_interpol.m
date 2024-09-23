function field_rescale = rescale_interpol(field,fx,gridsize,c)
% Utility function used while comparing differently sampled diffraction patterns (performs interpolation)
    % Dimension recalculation:
    span = -fx(1)+fx(end);  
    resolution = span / gridsize(1);
    start_index = round((c(1) + fx(end)) / resolution);
    end_index = round((c(end) + fx(end)) / resolution);

    field_cropped = field(start_index:end_index,start_index:end_index);

    [X_orig, Y_orig] = meshgrid(start_index:end_index, start_index:end_index);
    [X_new, Y_new] = meshgrid(linspace(start_index, end_index, length(c)), ...
                        linspace(start_index, end_index, length(c)));
    % length(c) == size of output matrix

    field_rescale = interp2(X_orig, Y_orig, field_cropped, X_new, Y_new, 'nearest');
end 
