#!/usr/bin/env python3
import sys
import math

def read_input():
    data = sys.stdin.read().strip().split()
    if not data:
        return 0, []
    it = iter(data)
    n = int(next(it))
    pts = []
    for _ in range(n):
        x = float(next(it)); y = float(next(it))
        pts.append((x,y))
    return n, pts

def polygon_area(pts):
    n = len(pts)
    a = 0.0
    for i in range(n):
        x1,y1 = pts[i]
        x2,y2 = pts[(i+1)%n]
        a += x1*y2 - x2*y1
    return abs(a) * 0.5

def polygon_perimeter_and_maxH(pts):
    n = len(pts)
    perim = 0.0
    maxH = float('inf')
    for i in range(n):
        x1,y1 = pts[i]
        x2,y2 = pts[(i+1)%n]
        L = math.hypot(x2-x1, y2-y1)
        perim += L
        # must keep L - 2H >= 0.1  => H <= (L - 0.1)/2
        maxH = min(maxH, (L - 0.1) / 2.0)
    if maxH < 0:
        maxH = 0.0
    return perim, maxH

def best_volume(pts):
    A = polygon_area(pts)
    P, maxH = polygon_perimeter_and_maxH(pts)
    best = 0.0
    # H must be a multiple of 0.1, starting at 0
    steps = int(math.floor(maxH * 10 + 1e-9))
    for k in range(steps + 1):
        H = k / 10.0
        inner = A - H * P + 4.0 * (H * H)
        if inner <= 0.0:
            continue
        vol = inner * H
        if vol > best:
            best = vol
    return best

def main():
    n, pts = read_input()
    if n == 0:
        print("0.00")
        return
    vol = best_volume(pts)
    print(f"{vol:.2f}")

if __name__ == "__main__":
    main()
