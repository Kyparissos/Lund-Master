function f = parzen(T,x)
    n = length(T);
    
    syms h
    for i = 1:n
        f = 0;
        x = (x-T(i))/h;
        t = exp(-x.^2./2)/(h*sqrt(2*pi));
        f = f+t;
    end
    f = f./n;
end
