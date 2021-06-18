function [] = save_multi_tiff(Im, name)

% Delete previous image
delete(name)

% Write multitiff
imwrite(Im(:,:,1)./max(max(Im(:,:,1))), name)
for i = 2:size(Im,3)
    imwrite(Im(:,:,i)./max(max(Im(:,:,i))), name, 'writemode', 'append')
end

