def HMS3PO(G,T,M,K,delta, eta, mu,epsilon,d,n):

    subsets = []
    for t in T:
        subset = {t}
        subsets.append(subset)
    
    left = G.keys() - set.union(*subsets)
    while left != set():
        shortest_dist = -1
        shortest_nodei = None
        shortest_nodeo = None
        shortest_set = None
        for unincorporated in left:
            l = [c[(unincorporated, dest)] for dest in Adj[unincorporated]]
            mini = argmin(l)
            dest = Adj[unincorporated][mini]
            if l[mini]<shortest_dist and dest in set.union(*subsets):
                shortest_dist = l[mini]
                shortest_nodei = unincorporated
                shortest_nodeo = Adj[unincorporated][mini]
                for subset in subsets:
                    if dest in subset:
                        shortest_set = subset
        if length(shortest_set) != d*len(G.V)/n-1:
            shortest_set.add(dest)
        else
            break
    
    

        

class G(V,E,conductivities):
    
    def init():
        self.V = V
        self.E = E
        self.conductivities = conductivities
        self.Adj = graph(V,E,c)
    
    def graph(V,E,c):
        Adj ={}

        for edge in E:
            origin = Adj.setdefault(edge[0],[])
            origin.append(edge[1])

            origin = Adj.setdefault(edge[1],[])
            origin.append(edge[0])

        return Adj
    
    def c(edge):
        return self.conductivities[edge]

def subset_graph(subset, destination, G):
    Adj ={}

    for key in G:
        if key in subset or key==destination:
        Adj[key] = G[key]

    return Adj
    

def bfs(G, start, end)
    visited = {} #no vertices visited thus far
    pathQueue = [[start]] #FIFO
    visited[start] = True
    while len(pathQueue) != 0:
        tmpPath = pathQueue.pop(0)
        lastNode = tmpPath[-1]
        if lastNode == end:
            return tmpPath
        for nextNode in graph.childrenOf(lastNode):
            if nextNode not in visited: #avoid visiting again
                newPath = tmpPath + [nextNode] #Note copy of list
                visited[nextNode] = True
                pathQueue.append(newPath)
    return None
