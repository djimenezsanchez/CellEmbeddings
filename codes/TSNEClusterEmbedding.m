function [Yembedding] = TSNEClusterEmbedding(Embedding,EmbeddingCluster,EmbeddingPerplexity,ColorComm)
addpath('FItSNE')
Yembedding = fast_tsne(Embedding);%,'Perplexity',EmbeddingPerplexity,'Exaggeration',12);         

figure;
gscatter(Yembedding(:,1),Yembedding(:,2),EmbeddingCluster,ColorComm); axis('tight');
ylabel('tsne dim. 1'); xlabel('tsne dim. 2')

