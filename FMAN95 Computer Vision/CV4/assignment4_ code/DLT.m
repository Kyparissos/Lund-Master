function X = DLT(x1, x2, P1, P2)
    X = [];
    for i=1:length(x1)
        M = [P1 -x1(:,i) zeros(3,1); P2 zeros(3,1) -x2(:,i)];
        [U,S,V] = svd(M);
        v = V(:,end);
        X=[X v(1:4,1)];
    end
    X = pflat(X);
end
