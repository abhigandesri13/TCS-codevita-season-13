def main():
    import sys
    from collections import deque, defaultdict

    data = sys.stdin.read().splitlines()
    n = int(data[0])
    segments = []
    graph = defaultdict(set)
    points_set = set()

    for i in range(1, n+1):
        x1, y1, x2, y2 = map(int, data[i].split())
        point1 = (x1, y1)
        point2 = (x2, y2)
        segments.append((point1, point2))
        points_set.add(point1)
        points_set.add(point2)
        graph[point1].add(point2)
        graph[point2].add(point1)

    visited = set()
    cycles = 0

    for point in points_set:
        if point not in visited and len(graph[point]) == 2:
            # Check if the component is a cycle
            queue = deque([point])
            comp = set()
            is_cycle = True
            while queue:
                current = queue.popleft()
                if current in visited:
                    continue
                visited.add(current)
                comp.add(current)
                if len(graph[current]) != 2:
                    is_cycle = False
                for neighbor in graph[current]:
                    if neighbor not in visited:
                        queue.append(neighbor)
            if is_cycle and len(comp) >= 3:
                cycles += 1
        # Also, if a node has degree 0 or 1, we mark it visited but skip.
        elif point not in visited:
            visited.add(point)

    print(cycles)

if __name__ == "__main__":
    main()