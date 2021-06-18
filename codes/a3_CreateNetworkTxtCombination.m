function [JsonEncodingG, JsonEncodingClassMap, JsonEncodingIdMap, AdjacencytoTxt, AttributestoTxt, maxNodeDegree,NeighborsToAnalysis] = a3_CreateNetworkTxt(SuperInfoCombination, ContextSize, mean_F,std_F)

AttrInfo = []; SuperInfoSize = []; SuperInfoColor = []; EdgeListAll = []; AdjacencytoTxt = []; NeighborsToAnalysis = [];
for i=1:size(SuperInfoCombination,2)
    % Extract Size and Color.
    SuperInfoSize = SuperInfoCombination{i}(:,3); 
    SuperInfoColor = [SuperInfoColor; SuperInfoCombination{i}(:,4:end)]; 
    
    % Find neighbors of Superpixels and put it in EdgeList
    AttrInfo = SuperInfoCombination{i}(:,1:2);
    
    
    Kneigh = 500; % numero de vecinos a encontrar 
    [cIdx,cD] = knnsearch(AttrInfo(:,1:2),AttrInfo(:,1:2),'K',Kneigh);
    cIdxSize = reshape(SuperInfoSize(cIdx),size(cIdx));
    CumcIdxSize = cumsum(cIdxSize,2); % Suma acumulada de los tamaños de los Superpixeles.
    MaskCumulative = CumcIdxSize<ContextSize; %Whether values are inside the context or not.
    MaskCumulative(:,1) = 1; % Have at least one edge. 
    cIdx(~MaskCumulative) = 0; %These values are eliminated later.
    NeighborsToAnalysis = [NeighborsToAnalysis; cIdx];
    cIdx = cIdx-1;
    % Max Node Degree.
    maxNodeDegree(i) = find(sum(MaskCumulative)==0,1, 'first'); %Actual
    % maxNodedegree
     [a,maxNodeDegree(i)] =min(diff(sum(MaskCumulative))); % Use of minimum of derivative.
    fprintf(['Number of max degrees:',num2str(maxNodeDegree(i)),'->',num2str(find(sum(MaskCumulative)==0,1, 'first'))]);

EdgeList = [];
for ii=2:find(sum(MaskCumulative)==0,1, 'first') %Avoid extra iterations 
    EdgeList = [EdgeList; [cIdx(:,1),cIdx(:,ii)]];         
end
% Edge to itself
EdgeList = [EdgeList; [[0:(size(SuperInfoCombination{i},1)-1)]', [0:(size(SuperInfoCombination{i},1)-1)]']];
% Eliminate Edges, if any of them is outside the context.
EdgeList = EdgeList(~sum(EdgeList==-1,2),:);
% Eliminate Edge Repetitions
EdgeList = sort(EdgeList,2);
EdgeList = unique(EdgeList,'rows');

% Send to Txt - Pytorch. In this case the first node is id=1.
AdjacencytoTxt = [AdjacencytoTxt; EdgeList+length(AdjacencytoTxt)+1]; % Adjacency matrix.

if i>1
    EdgeListAll = [EdgeListAll;EdgeList+max(EdgeListAll(:))+1];
else
    EdgeListAll = [EdgeList];
end

end
display(['Size of edge List: ',num2str(size(EdgeListAll))]);
maxNodeDegree = round(mean(maxNodeDegree));

AttributestoTxt = round((SuperInfoColor-mean_F)./std_F,5); % Node attributes.

% Encode JSON - Insert Graph attributes/edges in a struct.
% Save information for python
addpath('Numpy-Matlab');
% The header of the JSON, is it a directed graph?
Header = struct('directed',false,'graph',struct('name', 'disjoint_union( , )'),'multigraph',false);            
% Save Node information to JSON
columnHeadings = {'test' 'id' 'feature' 'val' 'label'};
% Preallocate structure    
nodes(1:size(SuperInfoColor,1)) = cell2struct(repmat({[]},numel(columnHeadings),1),columnHeadings,1);
% Add Node information iteratively
for ii=1:size(SuperInfoColor,1)              
        nodes(ii) = struct('test', false, 'id', num2str(ii-1), 'feature', round((SuperInfoColor(ii,:)-mean_F)./std_F,5), 'val',false, 'label',[0, 0]);
end
% Add links infor to the header
Header.nodes = nodes;
clear nodes;

% Save links information to JSON
columnHeadings = {'source' 'target' 'test_removed' 'train_removed'};
% Preallocate structure
links(1:size(EdgeListAll,1)) = cell2struct(repmat({[]},numel(columnHeadings),1),columnHeadings,1);
% Add Node information iteratively
for ii=1:size(EdgeListAll,1)
    links(ii) = struct('source', EdgeListAll(ii,1), 'target', EdgeListAll(ii,2), 'test_removed', false, 'train_removed', false);
end
% Add links information to the header
Header.links = links;
clear links;

% Generate encoding
JsonEncodingG = jsonencode(Header);

% Classmap with labels.
Id = cell(size(SuperInfoColor,1),1);
Label = cell(size(SuperInfoColor,1),1);
for ii=1:size(SuperInfoColor,1)
    Id{ii} = num2str(ii-1);
    Label{ii} = [ii-1];
end
JsonEncodingClassMap= jsonencode(containers.Map(Id, Label));

% IdMap with Ids of nodes. we can choose the name of the id.
Id = cell(size(SuperInfoColor,1),1);
NameId = cell(size(SuperInfoColor,1),1);
for ii=1:size(SuperInfoColor,1)
    Id{ii} = num2str(ii-1);
    NameId{ii} = ii-1;
end
JsonEncodingIdMap = jsonencode(containers.Map(Id, NameId));


end
