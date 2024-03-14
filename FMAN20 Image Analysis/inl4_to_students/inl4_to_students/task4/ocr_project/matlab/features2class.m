function Y = class_train(X,CD)
train_set = CD(1:size(CD,1)-1,:);% extract train data
diff = zeros(1,size(train_set,2));
% knn k=1
for i = 1:size(train_set,2)
    diff(i) = norm(X - train_set(:,i));% calculate the distance to each element
end
[value,position] = min(diff); % find index
Y = CD(end,position); % find label
end