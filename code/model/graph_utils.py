import numpy as np

def largest_scc(n, src, dst):
    """Iterative Kosaraju; returns a boolean mask of the largest strongly connected component."""
    out = [[] for _ in range(n)]
    inc = [[] for _ in range(n)]
    for a, b in zip(src.tolist(), dst.tolist()):
        out[a].append(b)
        inc[b].append(a)
    seen = [False] * n
    order = []
    for s in range(n):
        if seen[s]:
            continue
        seen[s] = True
        stack = [(s, 0)]
        while stack:
            x, i = stack.pop()
            if i < len(out[x]):
                stack.append((x, i + 1))
                y = out[x][i]
                if not seen[y]:
                    seen[y] = True
                    stack.append((y, 0))
            else:
                order.append(x)
    comp = [-1] * n
    sizes = []
    for s in reversed(order):
        if comp[s] >= 0:
            continue
        c = len(sizes)
        comp[s] = c
        stack = [s]
        cnt = 0
        while stack:
            x = stack.pop()
            cnt += 1
            for y in inc[x]:
                if comp[y] < 0:
                    comp[y] = c
                    stack.append(y)
        sizes.append(cnt)
    best = int(np.argmax(sizes))
    return np.array(comp) == best
