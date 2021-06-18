function [SuperpixelLabels, Superpixel_Colors, Superpixel_stats, numPixls] = a2_GenerateSuperpixel(Im,OriginalIm,Threshold)

% Load hierarchy superpixel algorithm
Compile_cpp_file = true; % Compile cpp file
if Compile_cpp_file; mex codes/HierarchySuperpixel/MultiSuperpixelHierarchyMex.cpp; end

% Image normalization
Im = Im-min(Im(:));
Im = uint16(round(Im/max(Im(:)).*(2^16-2)+1));    
OriginalIm = uint16(OriginalIm/max(OriginalIm(:)).*(2^16-2)+1);

% Hierarchy superpixel parameters
connectivity = 8; iter_switch = 0;
Perc_toJoin = 10; maxDistColor = 1000;
maxSize = 10000; doEuclidean = 1;
doMeanJoin = 1;

% Calculate hierarchy superpixel tree
tic;sh=MultiSuperpixelHierarchyMex(Im,Im,8,iter_switch,Perc_toJoin,maxDistColor,maxSize,doEuclidean,doMeanJoin);  
numPixls = double(sh.nvertex - sh.JoinedPixels);
disp(num2str(numPixls))
GetSuperpixels(sh,numPixls); [Superpixel_Colors,Superpixel_stats] = MeanColor(double(Im),double(sh.label)); 
SuperpixelLabels = sh.label+1;

% Find optimal number of superpixels
numPixls =round(numel(Im(:,:,1))*0.02);
for ii = 1:50
    GetSuperpixels(sh,numPixls); [Superpixel_Colors,Superpixel_stats] = MeanColor(double(Im),double(sh.label)); 
    Result = (sqrt(sum(sum(sum(abs(Superpixel_Colors-double(Im)).^2)))/numel(Im)))/(655.35*sqrt(size(Im,3)));%/(/(65535*sqrt(size(Im,3))));
    if Result>Threshold
        break
    end
    numPixls = ceil(numPixls/1.2);
end
SuperpixelLabels = sh.label+1;
display_multiplex(Superpixel_Colors,'Superpixel image')
fprintf(['\nNumber of superpixels: ',num2str(numPixls)]);
    


