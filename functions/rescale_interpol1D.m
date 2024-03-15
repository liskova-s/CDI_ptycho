function field_rescale = rescale_interpol1D(field,fx,gridsize,c)
    % Dimension recalculation:
    span = -fx(1)+fx(end);  
    resolution = span / gridsize(1);
    start_index = round((c(1) + fx(end)) / resolution);
    end_index = round((c(end) + fx(end)) / resolution);

    field_cropped = abs(field(start_index:end_index));

    X_orig = linspace(start_index,end_index, length(field_cropped));
    X_new = linspace(start_index, end_index, length(c));
    % length(c) == size of output matrix
    field_rescale = (interp1(X_orig,field_cropped, X_new,'nearest'))';
end 