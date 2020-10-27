# Note on Frontier and Path Decomposition

Consider a family of edge-induced subgraphs of a given graph.

## Path Decomposition

A path decomposition of a undirected graph $G=(V, E)$ is a path $\mathcal{P}$ satisfying the following conditions:
- (P1) every node of $\mathcal{P}$ is a subset of $V$.
- (P2) for any edge $\{ u, v \} \in E$, there exists a node $X \in \mathcal{P}$ such that $u, v \in X$.
- (P3) for any vertex $v \in V$, nodes of $\mathcal{P}$ containing $v$ induce sub-paths of $\mathcal{P}$.

By (P2), a path decomposition is derived when we fix the first appearance order of edges in $\mathcal{P}$.
```python
# construct a path decomposition and a corresponding frontiers with a fixed edge-order
from copy import deepcopy

path_decompositions = []
path_decompositions.append(set([]))
frontiers = []
bag = set([])
forget = set([])
for e in edges:
    # introduce vertices associated with edge `e`
    bag.insert(e.src)
    bag.insert(e.dst)
    # just after evaluating edge `e`
    frontiers.append(deepcopy(bag))
    path_decomposition.append(deepcopy(bag))
    # forget
    v_forget = []
    if neighbors(e.src).issubset(bag.union(forget)):
        v_forget.append(e.src)
    if neighbors(e.dst).issubset(bag.union(forget)):
        v_forget.append(e.dst)
    for v in v_forget:
        forget.insert(v)
        bag.remove(v)
    if v_forget:
        path_decomposition.append(deepcopy(bag))
```

A frontier is regarded as bag of a path decomposition just after introducing new vertexs.

```python
def get_child(state, edge, value):
    introduce(state, edge.src)
    introduce(state, edge.dst)
    process_edge(state, edge)
    forget(state)
```
