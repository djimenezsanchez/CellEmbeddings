function [Im] = ImageNormalization(Im)

    % Show original image
    display_multiplex(Im,'Original Image')
    
    % Calculate raw image statistics 
    ImageVec_Norm = reshape(Im,[size(Im,1)*size(Im,2),size(Im,3)]);
    ImageVec_NormMean = mean(ImageVec_Norm);
    ImageVec_NormStd = std(ImageVec_Norm);
    maxx = max(ImageVec_Norm);
    
    % Normalize image using z-score
    normalizeZSCORE=true;
    if normalizeZSCORE
        
        % Adjust image levels     
        ImageVec_Norm = reshape(Im,[size(Im,1)*size(Im,2),size(Im,3)]);    
        Im = reshape(ImageVec_Norm,size(Im));
        for i=1:size(Im,3);
            Im(:,:,i) = imadjust(Im(:,:,i)./maxx(i),stretchlim(Im(:,:,i)./maxx(i),[.001 .999]),[]);
            Im(:,:,i) = Im(:,:,i).*maxx(i);
        end    
        
        % Z-score normalization
        ImageVec_Norm = reshape(Im,[size(Im,1)*size(Im,2),size(Im,3)]);
        ImageVec_Norm = zscore(ImageVec_Norm);
       
%         % A sort of L2 normalization
%         ImageVec_Norm = (ImageVec_Norm'./sqrt(sum(ImageVec_Norm'.^2)))';
%         ImageVec_Norm(ImageVec_Norm<0) = -(sqrt(abs(ImageVec_Norm(ImageVec_Norm<0))+1)-1);  
%         ImageVec_Norm = ImageVec_Norm-min(ImageVec_Norm);

        figure; scatter_kde(ImageVec_Norm(1:100000:end,3),ImageVec_Norm(1:100000:end,4))
        figure; scatter_kde(ImageVec_Norm(1:100000:end,1),ImageVec_Norm(1:100000:end,2))     
        
    end
      
    % Show normalized image
    Im = reshape(ImageVec_Norm,[size(Im,1),size(Im,2),size(Im,3)]);
    display_multiplex(Im,'Normalized image')
end