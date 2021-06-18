clear all;
close all;
clc;
%% Unsupervised Learning of Contextual Information in Multiplex Immunofluorescence Tissue Cytometry
% Author: Daniel Jiménez-Sánchez
addpath(genpath('codes'));
LoadfolderName = 'Images/Original/';
SuperpixelfolderName = 'Images/Superpixel/';
Contextual_Cell_Emb_folderName = 'Images/Embedding/';
filenames = dir('Images/Original/*.TIFF');

%------- Options --------
HierarchySuperpixel_Thrs=2; % Threshold to select the number of total superpixels. 
                            % 3 is a good starting point. (Try with 2 or 4)      
ContextSize=150*150; % measured in pixels
%------------------------

%% Image Load and Normalization
Im=[];
for iIm=1:length(filenames)    
    % Load and Convert Image to double
    Imiter = double(loadTiff([LoadfolderName,filenames(iIm).name])); 
    Im = [Im,Imiter]; 
end 
% Image normalization
[Im] = a1_ImageNormalization(Im);
% Save Image.   
save_multi_tiff(Im,[SuperpixelfolderName,'Normalized_',filenames(iIm).name]);

%% Generate Superpixel Image.    
ImageVec_Norm = reshape(Im,[size(Im,1)*size(Im,2),size(Im,3)]);
ImageVec_Norm(ImageVec_Norm==0) = 10000*eps;       

% Hierarchy superpixel
addpath('Numpy-Matlab');
[SuperpixelLabels_All, Superpixel_Colors_All, Superpixel_stats_All,...
    numPixl] = a2_GenerateSuperpixel(Im, Im,HierarchySuperpixel_Thrs);

% Save Hierarchy superpixel labels
a2_save_HS_labels(LoadfolderName, SuperpixelfolderName, filenames, SuperpixelLabels_All);
clear SuperpixelLabels_All ImageSuper_All ImageSuperRGB_All SuperInfo_All ImageVec_Norm
    
%% Extract features 
npyFeatures = []; npySizes = [];
Superpixel_stats_all = {};
Superpixel_size_all = [];
for iIm=1:length(filenames)
    % Load Superpixel labels and original image.
    file_name = strsplit(filenames(iIm).name,'.'); file_name = strjoin(file_name(1:end-1),'.');
    load([SuperpixelfolderName,'SuperpixelLabels_',file_name,'.mat'],'SuperpixelLabels');
    Orig_Im = double(loadTiff([LoadfolderName,filenames(iIm).name]));
    
    % Extract color information. Superpixel_stats: [x_position, y_position, area, marker_mean_expression...]
    [~, Superpixel_stats] = MeanColor(Orig_Im, SuperpixelLabels);
    Superpixel_Morph = regionprops(SuperpixelLabels,'Eccentricity','ConvexArea','Circularity');
    Superpixel_stats = [Superpixel_stats,[Superpixel_Morph.Eccentricity]',...
                        [Superpixel_Morph.ConvexArea]',[Superpixel_Morph.Circularity]'];
    Superpixel_stats(Superpixel_stats==inf)=10;
    Superpixel_stats_all{iIm} = Superpixel_stats;
    Superpixel_size_all(iIm) = size(Superpixel_stats,1);
    npyFeatures = [npyFeatures; Superpixel_stats(:,4:end)];
    npySizes = [npySizes; Superpixel_stats(:,3)];
     
end    

%% Create graph of interconnected cells/superpixels
% Create Network with Nodes(Superpixels), Edges(Connections) and Attributes(Colors)    
[JsonEncodingG, JsonEncodingClassMap, JsonEncodingIdMap, AdjacencytoTxt, AttributestoTxt, maxNodeDegree, NeighborsToAnalysis] = a3_CreateNetworkTxtCombination(Superpixel_stats_all,ContextSize, mean(npyFeatures,1),std(npyFeatures,1));

% Save graph information 
mkdir([Contextual_Cell_Emb_folderName,'Graph_',num2str(ContextSize)]);    
csvwrite([Contextual_Cell_Emb_folderName,'Graph_',num2str(ContextSize),'/Superpixel_A.txt'],AdjacencytoTxt);    
csvwrite([Contextual_Cell_Emb_folderName,'Graph_',num2str(ContextSize),'/Superpixel_node_attributes.txt'],AttributestoTxt);    
fid = fopen([Contextual_Cell_Emb_folderName,'Graph_',num2str(ContextSize),'/G.json'], 'w');            
fwrite(fid, JsonEncodingG, 'char'); fclose(fid);
fid = fopen([Contextual_Cell_Emb_folderName,'Graph_',num2str(ContextSize),'/class_map.json'], 'w');            
fwrite(fid, JsonEncodingClassMap, 'char'); fclose(fid);        
fid = fopen([Contextual_Cell_Emb_folderName,'Graph_',num2str(ContextSize),'/id_map.json'], 'w');            
fwrite(fid, JsonEncodingIdMap, 'char'); fclose(fid);    
writeNPY(zscore(npyFeatures),[Contextual_Cell_Emb_folderName,'Graph_',num2str(ContextSize),'/feats.npy']);
writeNPY(npySizes,[Contextual_Cell_Emb_folderName,'Graph_',num2str(ContextSize),'/sizes.npy']);        

%% Encode cellular contextual information   
% Options to generate a embedding for each cell in the graph
model = 'graphsage_mean';
learning_rate = 0.010000; % Put 6 decimals to read correctly;
model_size = 'small';
max_degree = 60;
samples_1 = 60;
samples_2 = 30;
samples_3 = 30;
samples_4 = 30;
context_area = ContextSize;
dim_1 = 100;
dim_2 = 100;
hops = 2;
neg_sample_size = 30;
epochs =20;
batch_size=1500;
ParamFolder = [model,'_',model_size,'_',num2str(learning_rate,'%0.6f')];

% Run GraphSAGE    
commandStr = ['conda activate tf_gpu & ', ...
              'python -m graphsage.unsupervised_train --train_prefix Graph_',num2str(ContextSize), ...
              '/ --model ',model,' --learning_rate ',num2str(learning_rate),' --model_size ',model_size, ...
              ' --context_area ',num2str(context_area),' --max_degree ',num2str(max_degree),' --samples_1 ',num2str(samples_1),...
              ' --samples_2 ',num2str(samples_2),' --samples_3 ',num2str(samples_3),' --samples_4 ',num2str(samples_4),' --epochs ',num2str(epochs),...
              ' --neg_sample_size ',num2str(neg_sample_size),' --hops ',num2str(hops),' --batch_size ',num2str(batch_size),...
              ' --dim_1 ',num2str(dim_1), ' --dim_2 ',num2str(dim_2),...
              ' & pause'];

fileID = fopen('Images/Embedding/DoGraphSAGEembedding.bat','w');
fprintf(fileID,commandStr); fclose(fileID);

% Execute GraphSAGE
winopen([Contextual_Cell_Emb_folderName,'DoGraphSAGEembedding.bat']);
% delete([Contextual_Cell_Emb_folderName,'Graph_',num2str(ContextSize),'/',ParamFolder,'/val.npy'])
while not(isfile([Contextual_Cell_Emb_folderName,'Graph_',num2str(ContextSize),'/',ParamFolder,'/val.npy']))     
    pause(1)
end  

%% Cluster cell embeddings to infer spatial phenotypes.

% Load cell embeddings
H = readNPY([Contextual_Cell_Emb_folderName,'Graph_',num2str(ContextSize),'/',ParamFolder,'/val.npy']);
Hidx = load([Contextual_Cell_Emb_folderName,'Graph_',num2str(ContextSize),'/',ParamFolder,'/val.txt']);
clear Embeddingnow; Embedding([Hidx+1],:) = H;

% Cluster cell embeddings
n_clusters = 5;
CellEmbedding_Cluster = kmeans(Embedding(:,1:end),n_clusters,'MaxIter',200,'Replicates',10);

% Display Tsne of cell embeddings and inferred spatial phenotypes.
[Yembedding] = TSNEClusterEmbedding(zscore(Embedding')',CellEmbedding_Cluster(1:1:end,:),ceil(size(Embedding,1)*0.001),jet(4));

% Visualize spatial phenotypes on the images.
n_cell = 0;
for iIm=1:length(filenames)    
    % Load Superpixel label image
    file_name = strsplit(filenames(iIm).name,'.');
    file_name = strjoin(file_name(1:end-1),'.');    
    load([SuperpixelfolderName,'SuperpixelLabels_',file_name,'.mat']);
    
    % Visualize spatial phenotypes
    cec = CellEmbedding_Cluster(n_cell+1:n_cell+max(SuperpixelLabels(:)));
    SuperpixelLabel_clusters = cec(SuperpixelLabels);
    figure; imshow(SuperpixelLabel_clusters,[]); colormap(jet(4));
    save([SuperpixelfolderName,'Clustering_',file_name,'.mat'], 'SuperpixelLabel_clusters');
    imwrite(SuperpixelLabel_clusters./max(SuperpixelLabel_clusters(:)), [SuperpixelfolderName,'Clustering_',file_name,'.tif'])
    
end                     

