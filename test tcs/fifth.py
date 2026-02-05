def solve():
    import sys
    input = sys.stdin.read
    data = input().splitlines()

    N = int(data[0])
    lines = []
    for i in range(1, N + 1):
        x1, y1, x2, y2 = map(int, data[i].split())
        lines.append((x1, y1, x2, y2))
    K = int(data[N + 1])

    # ---------- Supercover Line ----------
    def supercover_line(x1, y1, x2, y2):
        """Return the set of grid cells touched by the line segment."""
        cells = set()

        dx = abs(x2 - x1)
        dy = abs(y2 - y1)
        x, y = x1, y1

        n = 1 + dx + dy
        x_inc = 1 if x2 > x1 else -1 if x2 < x1 else 0
        y_inc = 1 if y2 > y1 else -1 if y2 < y1 else 0
        error = dx - dy
        dx *= 2
        dy *= 2

        for _ in range(n):
            cells.add((x, y))
            if error > 0:
                x += x_inc
                error -= dy
            elif error < 0:
                y += y_inc
                error += dx
            else:  # exactly on a corner
                x += x_inc
                y += y_inc
                error += dx - dy
        return cells

    def cells_touched(x1, y1, x2, y2):
        return len(supercover_line(x1, y1, x2, y2))

    # ---------- Geometry helpers ----------
    def on_segment(x1, y1, x2, y2, x, y):
        """Check if (x,y) lies on the segment (x1,y1)-(x2,y2)."""
        return (min(x1, x2) <= x <= max(x1, x2)) and (min(y1, y2) <= y <= max(y1, y2))

    def intersect(l1, l2):
        """Find intersection point of two line segments (if any)."""
        x1, y1, x2, y2 = l1
        x3, y3, x4, y4 = l2

        denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
        if denom == 0:
            return None  # parallel

        px = ((x1*y2 - y1*x2)*(x3 - x4) - (x1 - x2)*(x3*y4 - y3*x4)) / denom
        py = ((x1*y2 - y1*x2)*(y3 - y4) - (y1 - y2)*(x3*y4 - y3*x4)) / denom

        if on_segment(x1, y1, x2, y2, px, py) and on_segment(x3, y3, x4, y4, px, py):
            return (round(px, 6), round(py, 6))
        return None

    # ---------- Step 1: Find all intersections ----------
    stars = {}
    for i in range(N):
        for j in range(i + 1, N):
            p = intersect(lines[i], lines[j])
            if p:
                stars.setdefault(p, set()).update([i, j])

    # ---------- Step 2: Compute star intensities ----------
    total_intensity = 0
    for p, line_ids in stars.items():
        if len(line_ids) == K:
            lengths = []
            px, py = p
            for lid in line_ids:
                x1, y1, x2, y2 = lines[lid]
                if (px, py) == (x1, y1) or (px, py) == (x2, y2):
                    # Case 1: star at endpoint → one segment
                    lengths.append(cells_touched(px, py, x2, y2))
                else:
                    # Case 2: star inside line → both halves
                    l1 = cells_touched(px, py, x1, y1)
                    l2 = cells_touched(px, py, x2, y2)
                    lengths.extend([l1, l2])
            intensity = min(lengths)
            total_intensity += intensity

    print(total_intensity)
