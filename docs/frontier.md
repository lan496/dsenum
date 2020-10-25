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
    # just before evaluating edge `e`
    frontiers.append(deepcopy(bag))
    # introduce vertices associated with edge `e`
    bag.insert(e.src)
    bag.insert(e.dst)
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

top-down construction
```python
# state_i: frontiers[i] -> Any
def get_child(i, state_i, value):
    e: Edge = edges[i]
    do_some_branching(state_i, e, value)

    init(state_next)
    for u in frontiers[i + 1]:
        update(state_next[u])

    # state_next: frontiers[i + 1] -> Any
    return state_next
```
