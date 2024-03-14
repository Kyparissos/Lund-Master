function dldX = relu_backward(X, dldY)
    dldX = (X>0).* dldY;
end
