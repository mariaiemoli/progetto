function [A] = ff2matlab(filename)
    fd = fopen(filename);
    % skip first three lines
    fgets(fd);
    fgets(fd);
    fgets(fd);
    % read the matrix
    [dim] = fscanf(fd, '%d %d %d %d', [1,4]);
    entries = fscanf(fd, '%f %f %f', [3, dim(4)]);
    A = sparse(entries(1,:), entries(2,:), entries(3,:), dim(1), dim(2), dim(4));
    A = full(A);
end