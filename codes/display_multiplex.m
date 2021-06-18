function display_multiplex(Im,name)

% Convert multiplex image to rgb
ImageRGB = zeros([size(Im(:,:,1)),3]); 
if size(Im,3)>=1
    ImageRGB(:,:,1) = ImageRGB(:,:,1) + Im(:,:,1); % First channel/marker to red
end
if size(Im,3)>=2
    ImageRGB(:,:,2) = ImageRGB(:,:,2) + Im(:,:,2); % Second channel/marker to green
end
if size(Im,3)>=3
    ImageRGB(:,:,3) = ImageRGB(:,:,3) + Im(:,:,3); % Third channel/marker to blue
end
if size(Im,3)>=4
    ImageRGB(:,:,1) = ImageRGB(:,:,1) + Im(:,:,4); 
    ImageRGB(:,:,2) = ImageRGB(:,:,2) + Im(:,:,4); % Fourth channel/marker to Yellow
end
if size(Im,3)>=5
    ImageRGB(:,:,1) = ImageRGB(:,:,1) + Im(:,:,5); 
    ImageRGB(:,:,3) = ImageRGB(:,:,3) + Im(:,:,5); % Fifth channel/marker to magenta
end
if size(Im,3)>=6
    ImageRGB(:,:,2) = ImageRGB(:,:,2) + Im(:,:,7); 
    ImageRGB(:,:,3) = ImageRGB(:,:,3) + Im(:,:,7); % sixth channel/marker to magenta
end  
    
% Normalize RGB colors to improve visualization 
for i=1:size(ImageRGB,3)
    ImageRGB(:,:,i)= ImageRGB(:,:,i)-min(min(ImageRGB(:,:,i))); 
    ImageRGB(:,:,i) = ImageRGB(:,:,i)./max(max(ImageRGB(:,:,i)));
end    
figure; imshow(ImageRGB,[]); title(name);    
