function f = my_parzen(x,T,h)
    n = size(T, 1);
    f = zeros(length(x), 1);
    for i = 1:length(x)
        xp = x(i);
        total = 0;
        for j = 1:n
            xbar = (xp - T(j))/h;
            K = exp(-xbar.^2/2)/sqrt(2*pi);
            total = total + K/h;
        end
        f(i) = total / n;
    end
end