clear;

Map_matrix = xlsread('ER_test1.xls');

% get the three centralities
[deg_cent, betw_cent, pagerank] = centrality_ranking(Map_matrix);

% get the steiner tree thing
[idx, set] = steiner_ranking(Map_matrix);

% get index of max for each centrality
[M_deg,I_deg] = max(deg_cent)
[M_betw,I_betw] = max(betw_cent)
[M_pgr,I_pgr] = max(pagerank)
% steiner tree performance vec
[M_st,I_st] = max(idx) 

p = plot(graph(set))
highlight(p, [I_deg, I_betw, I_pgr, I_st], 'NodeColor', 'r')
labelnode(p, I_deg, 'Highest Degree Centrality')
labelnode(p, I_betw, 'Highest Betweenness Centrality')
labelnode(p, I_pgr, 'Highest PageRanking and Degree Centrality')
labelnode(p, I_st, 'Best Steiner Tree')
title('Experimental Network with Highest Ranking Vertices Highlighted')
% highlight(p, I_betw)
% highlight(p, I_pgr)
% highlight(p, I_st)