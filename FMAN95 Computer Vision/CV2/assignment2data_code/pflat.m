function x = pflat(X);

[rows,cols] = size(X);
C_x = X((1:rows),:);
W = X(end, :);
x = C_x ./ W;


