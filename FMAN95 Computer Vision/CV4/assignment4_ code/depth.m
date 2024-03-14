function infront = depth(P1, P2, X)
    f1 = [];
    f2 = [];
    infront = 0;
    for i = 1:size(X,2)
        f1(:,i) = (sign(det(P1(:,1:3))))/(norm(P1(:,1:3),2)*X(4,i))*P1(3,:)*X(:,i);
        f2(:,i) = (sign(det(P2(:,1:3))))/(norm(P2(:,1:3),2)*X(4,i))*P2(3,:)*X(:,i);
        if sign(f1) > 0 & sign(f2) > 0 
            infront = infront + 1;
        end
    end
end
