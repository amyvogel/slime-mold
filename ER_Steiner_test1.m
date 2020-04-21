clear all; 
clc;

%{
Basically this file is my attempt to turn the Physarum inspired 
Multiple Source Single Sink algorithm for finding the Steiner Tree into
something that finds a Steiner Tree for every node in the set of terminals
with the condition that the tree is rooted at that node (I think idk I
could be way off with this code). What I did was take the part of the
original algorithm that selects the sink node randomly and changed it so
that the sink node is fixed and the algorithm optimizes a tree around that
given node. The idea is that if you choose some node to be the root of the
Steiner Tree and then compute the length of the tree based around it,
you'll be left with something that's like a centrality measure of how close
that one node is to all the other nodes in the set of terminals (or
foodpoints).

The evolutionary computation aspect of the code is wiped every
time a new node is chosen to be the root and there is no terminating
condition for optimality. The only thing that makes the algorithm move on
to testing the next node is the parameter kk, which is the number of
iterations.

In practice (at least with the network dataset given by Prof. 
Yahui Sun on github) several of the nine tested sink nodes reach the optimal 
Steiner Tree within 30 iterations, but my hope is that when testing more
nodes to be the roots, fewer root nodes will yield an optimal tree and the
idea of using this as a performance/centrality measure will make more sense
since the distribution of tree lengths will have a greater spread.

Honestly, I don't know how this code works and I mainly got this result by
fidgeting with the outer loops and deleting code. I really really hope that
some useable result can come out of this so the next step would be to test
it with a different network and see the distribution of tree lengths for
all the nodes.
%}

%{
The data file should be structured as:
    row 1: vertex number, NaN, NaN
    row 2: edge number, NaN, NaN

    row 3: first edge data (first node, second node, distance)
        ...
    row num_edges+2: last edge data (first node, second node, distance)

    row num_edges+3: NaN, NaN, NaN
    row num_edges+4: NaN, NaN, NaN
    row num_edges+5: NaN, NaN, NaN

    row num_edges+6: number of terminals, NaN, NaN

    row num_edges+7: first terminal node, NaN, NaN
        ...
    row num_edges+num_terminals+6: last terminal node, NaN, NaN
%}

%%%%%%%%%%%%%%%%%%%%%%%%% Physarum Optimization constructs a network    (one sink, many sources)
% filepath = 'D:\Python\Documents\MATLAB\ER_test1.csv';
filepath = 'ER_test1.csv';
opts=detectImportOptions(filepath);
data = readtable(filepath,opts);
data.Var1 = [];
data = data{:,:};

optimal=82;  %%% optimal length

fe_PO=0; % function evaluation
N=data(1,1); % vertex number
edge_num=data(2,1); 
terminal_num=data(edge_num+6,1);
foodpoint=zeros(N,1);
for i=(edge_num+7):(edge_num+terminal_num+6)
    terminal=data(i,1);
    foodpoint(terminal)=1;
end
foodnum=sum(foodpoint); 
%%%
food=zeros(foodnum,1);
r=0;
for i=1:N
      if foodpoint(i)==1
          r=r+1;
          food(r)=i;  %%%%%%  the num r food is the num i vertice
      end
end
%%%%%%%%%%%%%
for p=1:foodnum
    success=0; %%% success time
    solution=optimal+1;

    %%%%%%%%%%%%define flux
    I=1e0;  %%%%%
    kk=30;  %%% inner iteration times
    cutoff=1e-3;   %%%% cutoff value of D
    %%%%% evolution of conductivity
    alpha=0.122;
    sigma=1.0;
    %%%%%length matrix, set matrix, initial conductivity matrix
    v=1e4; %%% long distance
    L=v*ones(N,N);%%%
    set=zeros(N,N);%%%%%sets of adjacent points.   1 means it's adjacent, 0 means not
    D=zeros(N,N);
    edgevalue=ones(N,N);
    for i=3:(edge_num+2)
        L(data(i,1),data(i,2))=data(i,3); L(data(i,2),data(i,1))=data(i,3);
        set(data(i,1),data(i,2))=1; set(data(i,2),data(i,1))=1;
        D(data(i,1),data(i,2))=1; D(data(i,2),data(i,1))=1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Reduction tests
    % [Rset,RL,foodpoint,EdgeBefore,EdgeAfter,FoodB,FoodA,STedge,N,fuseL]=Reduction(solution,foodpoint,set,L);
    Rset=set; RL=L; fuseL=0;
    % foodnum=sum(foodpoint);
    % %%%
    % food=zeros(foodnum,1);
    % r=0;
    % for i=1:N
    %       if foodpoint(i)==1
    %           r=r+1;
    %           food(r)=i;  %%%%%%  the num r food is the num i vertice
    %       end
    % end
    %%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    whole=0;

    set=Rset; L=RL; D=set;
    %%%%%%%%%%%%%%%%%%%%%%%%whether a vertex exist
    exist=ones(N,1); %%%%%%%%%%%% 1 means exist, 0 means not

    %%% iteration of D and P
    for k=1:kk  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        NUM=zeros(N,1);
        %%%%%%%%%%%%%% whether a vertex exist
        linknum=zeros(N,1);  %% how many tubes linked to a certain vertex
        for i=1:N
            for j=1:N
                if set(i,j)==1
                    linknum(i)=linknum(i)+1;
                end
            end
            if linknum(i)==0
                exist(i)=0;
            end
        end
        vertex=sum(exist);  %%% number of existing vertex

        %%%%%%%%%%%%%%%%%%%%%%%%%%% unequal possibility 3
        adjacentlength=zeros(foodnum,1);
        for i=1:foodnum
             for j=1:N
                  if set(food(i),j)==1
                      adjacentlength(i)=adjacentlength(i)+L(food(i),j);
                  end
             end
        end
        [B,ind]=sort(adjacentlength);   %%%  ascending order,  B is the ordered vector, ind is the intex
        TotalAdjacent=0;
        for i=1:foodnum
            TotalAdjacent=TotalAdjacent+B(i);
        end

        sink = food(p);
        luck = p;

        %%%%%%%%%%%%%%   calculate the pressure
        A=zeros(vertex-1);  
        g=0; 
        for w1=1:N   
            if exist(w1)==1   %%% only calculate the existing vertice
                if w1==sink %%% neglect the vertex if it's the sink point
                    ;
                else
                    g=g+1; %% row number of A
                    h=0; 
                    NUM(w1)=g; %%%  record the row number's relationship with vertex number, which will be used when using pressure
                    for w3=1:N  %%%%%%%  the vertex corresponds to the columes of A
                        if exist(w3)==1  %%% only calculate the existing point
                             if w3==sink(1)  %%% neglect the vertex if it's the sink point
                                 ;
                             else
                                    h=h+1; %% colume number of A
                                    if w3==w1   %%%  the row vertex is the colume vertex
                                        for i=1:N
                                            if set(w1,i)==1
                                                A(g,h)=A(g,h)-D(w1,i)/L(w1,i);
                                            end
                                        end
                                    else    %%%%   the row vertex is not the colume vertex
                                        if set(w1,w3)==1
                                            A(g,h)=D(w1,w3)/L(w1,w3);
                                        end
                                    end
                             end
                        end
                    end
                end
            end
        end

        X=zeros(vertex-1,1); 
        %%%%%%%%%%%%% all foodpoints that are not sink are sources
        for i=1:foodnum
            if i==luck
                ;
            else
                X(NUM(food(i)))=-I;
            end
        end

        P=A\X;

        %%%%calculate the flux and update the conductivity
        for i=1:N
            for j=1:N
            if set(i,j)==1 
                if i==sink
                    Q(i,j)=D(i,j)/L(i,j)*(0-P(NUM(j)));
                        D(i,j)=edgevalue(i,j)*(D(i,j)+alpha*abs(Q(i,j))-sigma*D(i,j));
                        if D(i,j)<cutoff  
                          set(i,j)=0; set(j,i)=0;
                        end
                end
                if j==sink
                    Q(i,j)=D(i,j)/L(i,j)*(P(NUM(i))-0);
                        D(i,j)=edgevalue(i,j)*(D(i,j)+alpha*abs(Q(i,j))-sigma*D(i,j));
                        if D(i,j)<cutoff  
                          set(i,j)=0; set(j,i)=0;
                        end
                end
                if i~=sink && j~=sink
                    Q(i,j)=D(i,j)/L(i,j)*(P(NUM(i))-P(NUM(j)));
                        D(i,j)=edgevalue(i,j)*(D(i,j)+alpha*abs(Q(i,j))-sigma*D(i,j));
                        if D(i,j)<cutoff  
                          set(i,j)=0; set(j,i)=0;
                        end
                end
            end
            end
        end
        fe_PO=fe_PO+1; % function evaluation
        Total_length=0;  %%%  total length of the network
        for i=1:N
            for j=i:N
                if set(i,j)==1
                    Total_length=Total_length+L(i,j);
                end
            end
        end
        Total_length=Total_length+fuseL;
        fit_PO(fe_PO)=Total_length;

        EDGE(k)=sum(sum(set))/2;  %%%  record the evolution
        
        %%%% PLOTTING STEINER TREE SO FAR, ANIMATION %%%%%%
        % this is showing how the best steiner tree is found for each node
        figure(3);
        if k == 1
            baseSet = set;
            pl1 = plot(graph(baseSet));
        elseif k==kk
            pl1 = plot(graph(baseSet));
            highlight(pl1, graph(set));
            highlight(pl1, [p], 'NodeColor', 'r', 'EdgeColor', 'r', 'LineWidth', 3);
            pause(1);
        else
            pl1 = plot(graph(baseSet));
            highlight(pl1, graph(set));
            highlight(pl1, [p], 'NodeColor', 'r', 'EdgeColor', 'r');
            pause(.1);
        end
         
    end  %%%% k iteration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    min(fit_PO)

    %%%%%%%% plotting the tree %%%%%%%%
    figure(1);
    pl = plot(graph(Rset));
    highlight(pl, graph(set), 'LineWidth',3);
    highlight(pl, [p], 'NodeColor', 'r', 'MarkerSize', 6);
    pause(.33) % hold on a sec until next iteration, this is for the animation
    
    
    if Total_length==min(fit_PO)
        MinSet=set;
        for i=1:N
            for j=1:N
               if set(i,j)==1
                   edgevalue(i,j)=1.2;
               else 
                   edgevalue(i,j)=0.8;
               end
            end
        end
    end

    
end

%%%%%%%%%%%%%  Analysis
%%% note from Amy: I moved this out of the for-loop, didn't think it made a
%%% difference but it disrupted the animation.
figure(2);
plot(1:fe_PO,fit_PO)
grid on;

performance_vec = zeros(1,N);
j = 1;
for i = 1:fe_PO-1
    if fit_PO(i)<fit_PO(i+1)
        performance_vec(j) = fit_PO(j);
        j = j+1;
    end
end

performance_vec(50) = fit_PO(length(fit_PO));
performance_vec