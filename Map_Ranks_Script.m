%% Map 1 node ranking

Map_matrix = xlsread('ER_test1.xls');
[deg_cent, betw_cent, pagerank] = centrality_ranking(Map_matrix)

% now we will order the nodes from most important to least important
num_nodes = 50;
deg_order = [];
betw_order = [];
pagerank_order = [];
%degree centrality order
for node = 1:num_nodes
    % we compare each node to all other nodes already ordered 
    for position = 1:length(deg_order)
        if deg_cent(deg_order(position))>= deg_cent(node)
            break
        end
    end
    for position2 = 1:length(betw_order)
        if betw_cent(betw_order(position2))>= betw_cent(node)
            break
        end
    end
    for position3 = 1:length(pagerank_order)
        if pagerank(pagerank_order(position3))>= pagerank(node)
            break
        end
    end
    %we put the new node in the correct position in the ordering
    deg_order = [deg_order(1:position-1),node,deg_order(position:end)];
    betw_order = [betw_order(1:position2-1),node,betw_order(position2:end)];
    pagerank_order = [pagerank_order(1:position3-1),node,pagerank_order(position3:end)];
end

deg_order
betw_order
pagerank_order

        