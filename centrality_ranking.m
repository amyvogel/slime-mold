function [degree_centrality, betweeness_centrality , pagerank] = centrality_ranking(map)
% map must be a matrix, data must be read already
num_edges = map(2,2);
num_nodes = map(1,2);
endpoint = num_edges+2;

%get g raph characteristics
from = map(3:endpoint, 2)';
to = map(3:endpoint, 3)';
weight = map(3:endpoint, 4)';

%create graph
G = graph(from,to,weight);
pagerank = centrality(G, 'pagerank');
degree_centrality = centrality(G, 'degree');
betweeness_centrality = centrality(G, 'betweenness');
