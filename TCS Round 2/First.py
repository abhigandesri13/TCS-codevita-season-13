from collections import deque
import sys

def min_rooms_to_center(lines):
    N = len(lines)
    L = [len(s) for s in lines]
    open_rooms = set()
    for r, s in enumerate(lines):
        for i, ch in enumerate(s):
            if ch == '0':
                open_rooms.add((r, i))

    outer = N - 1
    dist = {}
    q = deque()

    for i in range(L[outer]):
        if (outer, i) in open_rooms:
            dist[(outer, i)] = 1
            q.append((outer, i))

    if not q:
        return -1

    while q:
        r, i = q.popleft()
        d = dist[(r, i)]

        if r == 0:
            return d

        # left & right
        for ni in ((i - 1) % L[r], (i + 1) % L[r]):
            if (r, ni) in open_rooms and (r, ni) not in dist:
                dist[(r, ni)] = d + 1
                q.append((r, ni))

        # inward
        inner_j = i // 2
        if (r - 1, inner_j) in open_rooms and (r - 1, inner_j) not in dist:
            dist[(r - 1, inner_j)] = d + 1
            q.append((r - 1, inner_j))

        # outward
        if r + 1 < N:
            for oi in (2*i, 2*i+1):
                oi %= L[r+1]
                if (r + 1, oi) in open_rooms and (r + 1, oi) not in dist:
                    dist[(r + 1, oi)] = d + 1
                    q.append((r + 1, oi))

    return -1

def main():
    data = sys.stdin.read().strip().split()
    N = int(data[0])
    lines = data[1:1+N]
    ans = min_rooms_to_center(lines)
    print(ans)   # <-- ONLY number, no newline issues

if __name__ == "__main__":
    main()
