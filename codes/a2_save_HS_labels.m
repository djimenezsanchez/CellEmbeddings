function [] = a2_save_HS_labels(LoadfolderName, SavefolderName, filenames, SuperpixelLabels_All)

i_image = 0;

% Save each superpixel label image in disk
for iIm=1:length(filenames)
    
    % Load original image
    ImOrig = double(loadTiff([LoadfolderName,filenames(iIm).name])); 
    
    % Extract superpixel labels
    SuperpixelLabels_original = SuperpixelLabels_All(:,i_image+1:i_image+size(ImOrig,1));  
    SuperpixelLabels = zeros(size(SuperpixelLabels_original));
    [uniq,a] = unique(SuperpixelLabels_original);
    itera = 1; uniqordered=[];
    for i=uniq'; uniqordered(i) = itera; itera=itera+1; end
    for x=1:size(SuperpixelLabels,1)
        for y=1:size(SuperpixelLabels,2)
            SuperpixelLabels(x,y) = uniqordered(SuperpixelLabels_original(x,y));
        end
    end
    
    % Save superpixel labels
    file_name = strsplit(filenames(iIm).name,'.');
    file_name = strjoin(file_name(1:end-1),'.');
    save([SavefolderName,'SuperpixelLabels_',file_name,'.mat'],'SuperpixelLabels');
    imwrite(SuperpixelLabels./max(SuperpixelLabels),[SavefolderName,'SuperpixelLabels_',filenames(iIm).name])
    
    i_image = i_image+size(ImOrig,1);
end
