function [F, STR] = get_features(image, mask)
%keyboard;
I  = image;
% Bi = imbinarize(I,60);
ed = edge(I,"canny");
SE = strel('disk',1);
ed = imdilate(ed,SE);

% F(3) = median(image(find(mask)));
% STR{3} = 'median intensity';

% F(1) = mean(image(find(mask)));
F(1) = mean(ed,[1 2]);
STR{1} = 'mean intensity';

% F(2) = std(image(find(mask)));
F(2) = std(ed,0,[1 2]);
STR{2} = 'std dev';
% 
% 
s = regionprops(ed,"Area","Extent","FilledArea","EulerNumber");
F(3) = s(1).Area;
STR{3} = 'Area';

F(4) = s(1).Extent;
STR{4} = 'Extent';
F(5) = s(1).FilledArea;
STR{5} = 'FilledArea';
F(6) = s(1).EulerNumber;
STR{6} = 'EulerNumber';



% F = extractHOGFeatures(ed);
% STR{1} = 'HOG';

% Need to name all features.
% assert(numel(F) == numel(STR));