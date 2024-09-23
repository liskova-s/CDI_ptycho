function c = generate_coordinates(gridsize, squaresize)
    % returns coordinate grid of gridsize elements with size squaresize
    % the origin (0,0) is situated in the middle of the grid (for odd
    % dimensions)
    rangex = (gridsize(1)-1)/2;
    x = (-rangex:rangex)*squaresize;
    rangey = (gridsize(1)-1)/2;
    y = (-rangey:rangey)*squaresize;
    [c(:,:,1),c(:,:,2)] = meshgrid(x, y);
end