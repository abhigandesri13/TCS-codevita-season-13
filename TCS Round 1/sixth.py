#!/usr/bin/env python3
import sys
import math
import itertools
from collections import deque, defaultdict

# ---------- Utilities to read input ----------
def read_ints():
    return list(map(int, sys.stdin.read().strip().split()))

# ---------- Build grid walls and find interior cells ----------
def build_walls(segments):
    # compute bounds
    xs = []
    ys = []
    for x1,y1,x2,y2 in segments:
        xs += [x1,x2]
        ys += [y1,y2]
    minx = min(xs); maxx = max(xs)
    miny = min(ys); maxy = max(ys)
    # we will create grid covering [minx, maxx] x [miny, maxy]
    W = maxx - minx
    H = maxy - miny
    # shift coords so min at 0
    shiftx = -minx
    shifty = -miny
    # walls: horizontal_walls[x][y] is wall between cell (x,y-1) and (x,y) at y line (0..H)
    horizontal = [[False]*(H+1) for _ in range(W)]
    # vertical_walls[x][y] is wall between cell (x-1,y) and (x,y) at x line (0..W)
    vertical = [[False]*H for _ in range(W+1)]
    for x1,y1,x2,y2 in segments:
        x1 += shiftx; x2 += shiftx; y1 += shifty; y2 += shifty
        if x1==x2:
            # vertical segment: from miny to maxy
            xa = x1
            ya = min(y1,y2)
            yb = max(y1,y2)
            for y in range(ya, yb):
                if 0 <= xa <= W and 0 <= y < H:
                    vertical[xa][y] = True
        elif y1==y2:
            ya = y1
            xa = min(x1,x2)
            xb = max(x1,x2)
            for x in range(xa, xb):
                if 0 <= x < W and 0 <= ya <= H:
                    horizontal[x][ya] = True
        else:
            # should not occur (inputs axis-aligned)
            pass
    return W, H, horizontal, vertical

def find_interior_cells(W, H, horizontal, vertical):
    # flood fill from outside: we consider bounding region expanded by 1 cell margin
    # grid cells indices 0..W-1, 0..H-1
    visited = [[False]*H for _ in range(W)]
    q = deque()
    # We'll BFS "outside" by considering positions outside bounding box and attempting to enter boundary cells where no wall
    # Start by pushing "virtual" neighbors around box edges: we will enqueue boundary cells which are reachable from outside across edges lacking walls.
    # enqueue any cell adjacent to an outer boundary where there is no wall on that outer boundary side.
    # top boundary y=H: horizontal[x][H] is wall at top edge. If no wall, then cell (x,H-1) is reachable from outside.
    for x in range(W):
        if H>0 and not horizontal[x][H]:
            if not visited[x][H-1]:
                visited[x][H-1] = True
                q.append((x,H-1))
        if H>0 and not horizontal[x][0]:
            if not visited[x][0]:
                visited[x][0] = True
                q.append((x,0))
    for y in range(H):
        if W>0 and not vertical[0][y]:
            if not visited[0][y]:
                visited[0][y] = True
                q.append((0,y))
        if W>0 and not vertical[W][y]:
            if not visited[W-1][y]:
                visited[W-1][y] = True
                q.append((W-1,y))
    # Also if grid empty, nothing to do
    # BFS moves: left,right,up,down if no wall between cells
    while q:
        x,y = q.popleft()
        # left: to x-1
        if x-1 >= 0:
            # wall between (x-1,y) and (x,y) is vertical[x][y]
            if not vertical[x][y] and not visited[x-1][y]:
                visited[x-1][y] = True
                q.append((x-1,y))
        # right: to x+1
        if x+1 < W:
            if not vertical[x+1][y] and not visited[x+1][y]:
                visited[x+1][y] = True
                q.append((x+1,y))
        # down: y-1
        if y-1 >= 0:
            if not horizontal[x][y] and not visited[x][y-1]:
                visited[x][y-1] = True
                q.append((x,y-1))
        # up: y+1
        if y+1 < H:
            if not horizontal[x][y+1] and not visited[x][y+1]:
                visited[x][y+1] = True
                q.append((x,y+1))
    # interior cells are those not visited but within 0..W-1,0..H-1 and that are inside polygon (i.e., enclosed)
    interior = []
    for x in range(W):
        for y in range(H):
            if not visited[x][y]:
                interior.append((x,y))
    return interior, W, H

# ---------- Generate triples A,B,C ----------
def find_candidate_triples(area):
    triples = set()
    # iterate A,B,C from 1..area
    # We can enumerate factors: but brute small
    for A in range(1, area+1):
        for B in range(1, area+1):
            for C in range(1, area+1):
                if 2*(A*B + A*C + B*C) == area:
                    triples.add(tuple(sorted((A,B,C))))
    return sorted(triples)

# ---------- Tiling interior by rectangles ----------
def cells_to_index_map(cells):
    idx = {c:i for i,c in enumerate(cells)}
    return idx

def neighbors4(cell):
    x,y = cell
    return [(x-1,y),(x+1,y),(x,y-1),(x,y+1)]

def tile_with_rectangles(cells_set, W, H, dims_multiset):
    # dims_multiset: list of (w,h) possibilities but we accept rectangles by area and allow orientation
    # We'll backtrack placing rectangles. Represent sheet coords by set cells_set.
    cells = sorted(cells_set)
    idx = cells_to_index_map(cells)
    occupied = {c: -1 for c in cells}  # -1 unassigned, else face id
    target_counts = dims_multiset[:]  # list of (area, w, h) to place counts
    total = len(cells)
    placements = []

    # prepare list of rectangles by area counting duplicates required (we will pick specific w,h possibilities per area)
    # dims_multiset is list of areas with multiplicities but we'll store allowed (w,h) for each area when trying placements

    # To speed up, we will create a list of areas to place (6 items)
    areas_list = [d for d in dims_multiset]  # each element is list of (w,h) possibilities for that face
    # But easier approach: we will create list of target areas values e.g. [AB,AB,AC,AC,BC,BC] and choose (w,h) fitting each later.
    # For clarity dims_multiset currently will be list of tuples (area,). We'll pass a mapping area->list of possible (w,h) dims.

def generate_rectangle_orientations_for_area(area, A,B,C):
    # Given area and A,B,C values, produce all (w,h) integer pairs that multiply to area and are among {A,B,C} combos
    # Valid rectangle side lengths must be two of (A,B,C)
    vals = [A,B,C]
    res = []
    for x,y in itertools.permutations(vals,2):
        if x*y == area:
            res.append((x,y))
    # remove duplicates
    res = list(dict.fromkeys(res))
    return res

# We'll implement a tiler that tries to tile using specific required (w,h) rectangles (6 of them).
def try_tile(cells_set, W, H, rect_sizes):
    # rect_sizes: list of (w,h) for 6 rectangles (order matters). We'll try all assignments/permutations outside if needed.
    cells = sorted(cells_set)
    cellset = set(cells)
    occupied = {}
    placements = [None]*6

    # Precompute bounding: we use coordinates as small ints (x,y)
    cells_sorted = sorted(cells, key=lambda c:(c[1],c[0]))  # by y then x to get top-left scanning

    def find_next_unassigned():
        for c in cells_sorted:
            if c not in occupied:
                return c
        return None

    def can_place_at(tid, x0,y0, w,h):
        # try to place rectangle with top-left at (x0,y0) spanning width w in +x and height h in +y (we'll allow placing in any orientation by choosing x0,y0 appropriately).
        # But our grid uses cell coordinates increasing x to right and y upward. All good.
        for dx in range(w):
            for dy in range(h):
                c = (x0+dx, y0+dy)
                if c not in cellset or c in occupied:
                    return False
        return True

    def place_rect(tid, x0,y0,w,h):
        coords = []
        for dx in range(w):
            for dy in range(h):
                c = (x0+dx, y0+dy)
                occupied[c] = tid
                coords.append(c)
        placements[tid] = coords

    def remove_rect(tid):
        coords = placements[tid]
        for c in coords:
            del occupied[c]
        placements[tid] = None

    # For each tid, we'll find candidate top-left positions by scanning all cells and attempting to fit w,h
    # But top-left meaning minimum x and minimum y among placed cells; we must ensure all cells fit inside grid
    # We'll find the first unassigned cell and try to place any remaining rectangle so that it covers that cell.
    used = [False]*6

    def backtrack(placed):
        if placed == 6:
            return True
        nxt = find_next_unassigned()
        if nxt is None:
            return False
        x0_cell,y0_cell = nxt
        # try any unused rectangle that can cover this cell by choosing a top-left x,y such that x0_cell in [x0, x0+w-1] and y0_cell in [y0, y0+h-1]
        for tid in range(6):
            if used[tid]: continue
            w,h = rect_sizes[tid]
            # try all possible top-left x0,y0 such that x0_cell in range(x0,x0+w) and y0_cell in range(y0,y0+h)
            for x0 in range(x0_cell - (w-1), x0_cell+1):
                for y0 in range(y0_cell - (h-1), y0_cell+1):
                    # quick bounds test: ensure rectangle within grid box
                    # we don't know global grid max but cells_set used; check can_place_at
                    if can_place_at(tid, x0, y0, w, h):
                        used[tid] = True
                        place_rect(tid,x0,y0,w,h)
                        if backtrack(placed+1):
                            return True
                        remove_rect(tid)
                        used[tid] = False
        return False

    ok = backtrack(0)
    if not ok:
        return None
    # Return placements mapping tid->list of cells
    return placements

# ---------- Build adjacency between rectangles ----------
def build_adjacency(placements):
    # placements: list of lists of cells for each face id
    face_count = len(placements)
    face_cells = [set(p) for p in placements]
    adj = [[] for _ in range(face_count)]
    # Two faces adjacent if they share an edge segment (i.e., there exists a cell in face A and neighbor cell in face B sharing side)
    for i in range(face_count):
        for j in range(i+1, face_count):
            found = False
            for (x,y) in face_cells[i]:
                for nx,ny in [(x-1,y),(x+1,y),(x,y-1),(x,y+1)]:
                    if (nx,ny) in face_cells[j]:
                        # they share an edge
                        adj[i].append(j)
                        adj[j].append(i)
                        found = True
                        break
                if found: break
    return adj

# ---------- Folding simulation ----------
# We'll represent faces as rigid rectangles with integer size in grid units and their polygon coordinates in 2D.
# We will map each rectangle to a 3D face on a cuboid with dims (A,B,C) by assigning for base face an orientation (u_axis -> one cuboid axis, v_axis -> another) and position.
# Each face placement will map its unit cells to unit square positions on cuboid surface identified by (face_id_on_cuboid, u, v) where face_id_on_cuboid is one of 6 faces.
# We'll perform BFS/DFS across adjacency and fold neighbors by rotating 90 degrees around shared edge.

# Cuboid faces indexing:
# 0: x=0 plane (vary y in [0,B), z in [0,C))  -- dims B x C
# 1: x=A plane (vary y,z)                     -- dims B x C
# 2: y=0 plane (vary x,z)                     -- dims A x C
# 3: y=B plane                                -- dims A x C
# 4: z=0 plane (vary x,y)                     -- dims A x B
# 5: z=C plane                                -- dims A x B

# For a placed face we store mapping for every 2D sheet cell -> (cub_face, u, v) integers.

def face_dims_for_cuboid(A,B,C):
    # returns dict face_id -> (w,h) dims in units for that face (u_range,h_range)
    return {
        0: (B, C),
        1: (B, C),
        2: (A, C),
        3: (A, C),
        4: (A, B),
        5: (A, B),
    }

# For mapping orientation, we will place base onto any of the 6 cuboid faces.
# For base mapping we need to choose which cuboid face and an orientation (flip/rotate) so that the rectangle dims match the face dims.

def try_fold(placements, adj, A,B,C):
    # placements: list of lists of (x,y) coords (sheet coordinate grid). We will compute each face's local coordinates relative to its rectangle's top-left
    face_count = len(placements)
    cellsets = [set(p) for p in placements]
    # Precompute bounding rectangles and local coordinate mapping (local origin chosen as min x and min y of face cells)
    face_bounds = []
    face_local = []
    for p in placements:
        xs = [c[0] for c in p]; ys = [c[1] for c in p]
        minx=min(xs); miny=min(ys)
        maxx=max(xs); maxy=max(ys)
        w = maxx-minx+1; h = maxy-miny+1
        face_bounds.append((minx,miny,w,h))
        # create set of local coords present
        local = set((x-minx,y-miny) for (x,y) in p)
        face_local.append(local)

    # For a candidate mapping, we will assign each face to a cuboid face index and place each face's local coordinates to (u,v) on that cuboid face.
    face_dims = face_dims_for_cuboid(A,B,C)

    # helper: check rectangle dims match target dims (either orientation)
    def dims_match(w,h, target):
        tw,th = target
        return (w==tw and h==th) or (w==th and h==tw)

    # create adjacency edges list
    # We'll choose any face as base and try mapping to any cuboid face and any orientation (rotations/reflections).
    # For orientation, for mapping local coords (u_local,v_local) to cuboid face coords (u,v), we can choose any of 4 transforms: (u,v),(v,tw-1-u),(tw-1-u,th-1-v),(th-1-v,u) - rotations and flips.
    def get_orientations(w,h,tw,th):
        # returns list of functions map_local(u,v)->(u2,v2) on target face dims tw x th
        # We'll return list of possible transforms as lambdas; but for speed return small set of explicit transforms.
        res = []
        # If w==tw and h==th then rotations allowed
        if w==tw and h==th:
            res.append(lambda u,v: (u,v))
            res.append(lambda u,v: (w-1-u, h-1-v))
            res.append(lambda u,v: (v, w-1-u))  # rotate 90 (works only if w==h?) but safe if w==h==tw==th
            res.append(lambda u,v: (h-1-v, u))
        elif w==th and h==tw:
            # swapped orientation
            res.append(lambda u,v: (v,w-1-u))
            res.append(lambda u,v: (w-1-v, u))
            res.append(lambda u,v: (u,v))
            res.append(lambda u,v: (w-1-u, h-1-v))
        else:
            # fallback: only direct placement if sizes equal
            res.append(lambda u,v:(u,v))
        return res

    # map to track occupied cuboid face unit squares -> which sheet face cell occupies it
    # Try base face assignments
    for base_face in range(face_count):
        w0,h0 = face_bounds[base_face][2], face_bounds[base_face][3]
        # for each target cuboid face id that has matching dims
        for cub_face in range(6):
            tw,th = face_dims[cub_face]
            if not dims_match(w0,h0,(tw,th)): continue
            # try possible orientations mapping local -> cuboid coords
            orients = get_orientations(w0,h0,tw,th)
            # But need to ensure orientation functions are valid; we will iterate them and check placement
            for orient in orients:
                # structure to hold assignments
                assigned = [None]*face_count  # for each face store (cub_face_id, offset_u, offset_v, orientation_func, tw,th)
                # occupancy: dict (cub_face,u,v) -> (faceid, localcoord)
                occ = {}
                # we will perform DFS spreading from base
                # For base, set offset_u,offset_v such that its top-left local (0,0) maps to some u0,v0 on cub_face so that all mapped coords within 0..tw-1,0..th-1
                # Try all possible offsets where base's mapped coordinates fit
                possible_offsets = []
                # base local max coords
                minx,miny,w,h = face_bounds[base_face]
                # compute mapped bounding for local coords (0..w-1,0..h-1)
                # For offset candidates u0 in [0..tw-w], v0 in [0..th-h] if mapping preserves orientation without swap.
                for u0 in range(0, tw - w + 1):
                    for v0 in range(0, th - h + 1):
                        # verify mapping within bounds
                        ok_flag = True
                        for (u,v) in face_local[base_face]:
                            uu,vv = orient(u,v)
                            U = u0 + uu; V = v0 + vv
                            if not (0 <= U < tw and 0 <= V < th):
                                ok_flag = False; break
                        if ok_flag:
                            possible_offsets.append((u0,v0))
                if not possible_offsets:
                    continue
                # now try each possible offset
                for u0,v0 in possible_offsets:
                    assigned = [None]*face_count
                    occ = {}
                    # place base
                    assigned[base_face] = (cub_face, u0, v0, orient, tw, th)
                    ok_place = True
                    for (u,v) in face_local[base_face]:
                        U,V = orient(u,v)
                        key = (cub_face, u0+U, v0+V)
                        if key in occ:
                            ok_place = False; break
                        occ[key] = (base_face, (u,v))
                    if not ok_place: continue
                    # BFS fold neighbors: maintain stack of placed faces
                    stack = [base_face]
                    success = True
                    while stack and success:
                        f = stack.pop()
                        f_cub, f_u0, f_v0, f_orient, f_tw, f_th = assigned[f]
                        # examine neighbors on sheet
                        for nb in adj[f]:
                            if assigned[nb] is not None: continue
                            # Determine shared edge(s) between f and nb in sheet local coords
                            # find any adjacent cell pair (a in f, b in nb) that are neighbors -> that edge is the fold line
                            shared_edges = []
                            for (x,y) in placements[f]:
                                for dx,dy,dir_flag in [(-1,0,'L'),(1,0,'R'),(0,-1,'D'),(0,1,'U')]:
                                    nx,ny = x+dx,y+dy
                                    if (nx,ny) in set(placements[nb]):
                                        # on f local coords (uF,vF) and nb local coords (uN,vN)
                                        minxf,minyf,fw,fh = face_bounds[f]
                                        minxn,minyn,wn,hn = face_bounds[nb]
                                        uF = x - minxf; vF = y - minyf
                                        uN = nx - minxn; vN = ny - minyn
                                        shared_edges.append(( (uF,vF),(uN,vN), dir_flag ))
                            if not shared_edges:
                                success = False; break
                            # pick any shared edge (there may be many along full edge); compute folding transform around that edge
                            # We'll take first shared edge
                            (uF,vF),(uN,vN),dir_flag = shared_edges[0]
                            # location on placed face in cuboid coords:
                            UF,VF = f_orient(uF,vF); # local mapping
                            cub_u = f_u0 + UF; cub_v = f_v0 + VF
                            # The shared edge in cuboid coords lies along one of the local face edges: determine edge direction vector in local face coords
                            # For f, the neighbor cell is at dir_flag, which indicates which side of f is the shared edge.
                            # Need to compute the target cub_face for nb after 90-degree fold.
                            # We deduce the fold by using the fact that both cells are adjacent in sheet; after folding, nb must occupy the face adjacent to f along that edge on the cuboid.
                            # To determine the neighbor cub_face, we find which cub_face side of f has that edge.
                            # For face f at cub_face id, edges correspond to directions along its u/v axes: we'll map edges to world: each cub_face has 4 possible neighboring faces.
                            # We'll create a small map of neighbor relationships for cuboid faces as function of edge (left,right,up,down).
                            def get_face_neighbours(cface, A,B,C):
                                # returns mapping {'L':(face_id, transform), ...} where transform maps coords from neighbor local to current local after folding? We'll instead produce adjacency topology only.
                                # We'll predefine neighbors for cuboid faces by inspection.
                                # face dims: 0 x=0 (B x C) u->y, v->z
                                # It's easier to hardcode neighbors: neighbors[cface][edge] -> (neighbor_face, how the neighbor face u/v coordinates correspond when unfolded)
                                # We'll instead compute neighbor faces via a static cube net reference.
                                # Hardcode using conventional orientation:
                                pass
                            # Rather than building complex neighbor mapping, we'll simulate geometry using vectors: represent each cub_face with origin point and two axis vectors in 3D; then folding is a rotation of neighbor plane around the shared edge (3D vector).
                            # Build canonical cuboid coordinate axes in 3D: x in [0,A], y in [0,B], z in [0,C]
                            # For each cub_face id, define its origin (corner) and local axes (u_dir,v_dir) mapping to 3D unit vectors and ranges.
                            def cuboid_face_frame(A,B,C, face_id):
                                # returns origin (ox,oy,oz), u_vec (dx,dy,dz), v_vec (dx,dy,dz), u_len, v_len
                                if face_id==0:
                                    # x=0, u=y(0..B-1), v=z(0..C-1)
                                    return (0,0,0), (0,1,0), (0,0,1), B, C
                                if face_id==1:
                                    # x=A, u=y increasing, v=z increasing, origin at x=A,y=0,z=0
                                    return (A,0,0), (0,1,0), (0,0,1), B, C
                                if face_id==2:
                                    # y=0, u=x increasing, v=z increasing
                                    return (0,0,0), (1,0,0), (0,0,1), A, C
                                if face_id==3:
                                    # y=B
                                    return (0,B,0), (1,0,0), (0,0,1), A, C
                                if face_id==4:
                                    # z=0, u=x increasing, v=y increasing
                                    return (0,0,0), (1,0,0), (0,1,0), A, B
                                if face_id==5:
                                    # z=C
                                    return (0,0,C), (1,0,0), (0,1,0), A, B
                            # Now compute the 3D coordinates of the shared edge on face f
                            f_origin, f_uvec, f_vvec, f_ulen, f_vlen = cuboid_face_frame(A,B,C,f_cub)
                            # compute 3D pos of the unit square at (UF,VF) -> it's the lower-left corner at (UF, VF) in local axes
                            fx = f_origin[0] + f_uvec[0]*UF + f_vvec[0]*VF
                            fy = f_origin[1] + f_uvec[1]*UF + f_vvec[1]*VF
                            fz = f_origin[2] + f_uvec[2]*UF + f_vvec[2]*VF
                            # the shared edge lies between this cell and neighbor cell: direction vector in local of f indicates which edge orientation
                            if dir_flag=='L':
                                edge_start = (fx,fy,fz)
                                edge_dir = (-f_uvec[0], -f_uvec[1], -f_uvec[2])
                            elif dir_flag=='R':
                                edge_start = (fx+f_uvec[0], fy+f_uvec[1], fz+f_uvec[2])
                                edge_dir = (f_uvec[0], f_uvec[1], f_uvec[2])
                            elif dir_flag=='D':
                                edge_start = (fx,fy,fz)
                                edge_dir = (-f_vvec[0], -f_vvec[1], -f_vvec[2])
                            else: # 'U'
                                edge_start = (fx+f_vvec[0], fy+f_vvec[1], fz+f_vvec[2])
                                edge_dir = (f_vvec[0], f_vvec[1], f_vvec[2])
                            # The neighbor face must be the face adjacent to f across that edge in the cuboid geometry.
                            # Find which cub_face has cells adjacent to these coordinates across that edge.
                            # We'll inspect all six faces and check if there's a face whose frame includes the adjacent cell across edge_dir direction by 0.5 unit.
                            neighbor_found = False
                            for candidate in range(6):
                                if candidate==f_cub: continue
                                c_origin,c_uvec,c_vvec,c_ulen,c_vlen = cuboid_face_frame(A,B,C,candidate)
                                # We need to check if there exists integer Uc,Vc such that the cell on candidate shares that edge.
                                # We'll try mapping nb local coords to 3D and see if any of them equals edge cell across direction.
                                # brute force: for each local cell in nb's local set, compute its 3D lower-left corner and see if it matches edge_start + small shift in direction.
                                foundmatch = False
                                for (u_nb,v_nb) in face_local[nb]:
                                    ux = c_origin[0] + c_uvec[0]*u_nb + c_vvec[0]*v_nb
                                    uy = c_origin[1] + c_uvec[1]*u_nb + c_vvec[1]*v_nb
                                    uz = c_origin[2] + c_uvec[2]*u_nb + c_vvec[2]*v_nb
                                    # If this cell is located adjacent to (fx,fy,fz) along edge_dir direction, then the vector from (ux,uy,uz) to (fx,fy,fz) should be +/- unit vec from uvec/vvec etc.
                                    # compute delta
                                    dx = fx - ux; dy = fy - uy; dz = fz - uz
                                    # check if this delta is one of unit directions (±uvec or ±vvec)
                                    okdelta = False
                                    for vec in [c_uvec, c_vvec]:
                                        if (abs(dx - vec[0])<1e-9 and abs(dy - vec[1])<1e-9 and abs(dz - vec[2])<1e-9) or (abs(dx + vec[0])<1e-9 and abs(dy + vec[1])<1e-9 and abs(dz + vec[2])<1e-9):
                                            okdelta = True
                                            break
                                    if okdelta:
                                        foundmatch = True
                                        break
                                if foundmatch:
                                    # candidate is neighbor face
                                    nb_face_id = candidate
                                    neighbor_found = True
                                    break
                            if not neighbor_found:
                                success = False
                                break
                            # Now compute placement for neighbor face: we need to find offset and orientation that map nb local coords to candidate face coords so that the shared edge matches and no overlaps
                            # We'll attempt all possible orientations and offsets on candidate face that fit nb's local coords into its face dims.
                            c_tw,c_th = face_dims[nb_face_id]
                            # compute orientation funcs for nb's local dims (wn,hn)
                            wn,hn = face_bounds[nb][2], face_bounds[nb][3]
                            orient_candidates = get_orientations(wn,hn,c_tw,c_th)
                            placed_nb = False
                            for or_nb in orient_candidates:
                                # try offsets u0,v0 within [0..c_tw-wn]x[0..c_th-hn]
                                for u0c in range(0, c_tw - wn + 1):
                                    for v0c in range(0, c_th - hn + 1):
                                        # check alignment of shared edge: for the shared local pair (uN,vN) we had earlier, its mapped 3D coords should be adjacent to f's mapped coords
                                        Uc,Vc = or_nb(uN,vN)
                                        # compute mapped 3D coords of nb cell lower-left
                                        c_origin2, c_uvec2, c_vvec2, _, _ = cuboid_face_frame(A,B,C,nb_face_id)
                                        nb3x = c_origin2[0] + c_uvec2[0]*(u0c+Uc) + c_vvec2[0]*(v0c+Vc)
                                        nb3y = c_origin2[1] + c_uvec2[1]*(u0c+Uc) + c_vvec2[1]*(v0c+Vc)
                                        nb3z = c_origin2[2] + c_uvec2[2]*(u0c+Uc) + c_vvec2[2]*(v0c+Vc)
                                        # difference between nb3 and f3 should be +/- one unit in direction of edge
                                        dx = fx - nb3x; dy = fy - nb3y; dz = fz - nb3z
                                        mag = abs(dx) + abs(dy) + abs(dz)
                                        if abs(mag - 1.0) > 1e-9:
                                            continue
                                        # Now check that all nb local cells map into candidate face coords and do not collide with occ
                                        collide = False
                                        for (u_loc,v_loc) in face_local[nb]:
                                            Upl,Vpl = or_nb(u_loc,v_loc)
                                            key = (nb_face_id, u0c+Upl, v0c+Vpl)
                                            if key in occ:
                                                collide = True; break
                                            # bounds check
                                            if not (0 <= u0c+Upl < c_tw and 0 <= v0c+Vpl < c_th):
                                                collide = True; break
                                        if collide:
                                            continue
                                        # accept placement
                                        assigned[nb] = (nb_face_id, u0c, v0c, or_nb, c_tw, c_th)
                                        for (u_loc,v_loc) in face_local[nb]:
                                            Upl,Vpl = or_nb(u_loc,v_loc)
                                            occ[(nb_face_id, u0c+Upl, v0c+Vpl)] = (nb, (u_loc,v_loc))
                                        stack.append(nb)
                                        placed_nb = True
                                        break
                                    if placed_nb: break
                                if placed_nb: break
                            if not placed_nb:
                                success = False
                                break
                        # end for neighbors
                    # end while stack
                    if success:
                        # Check occ covers exactly all cuboid surface unit squares (2*(AB+AC+BC) units)
                        expected = 2*(A*B + A*C + B*C)
                        if len(occ) == expected:
                            # ensure occ mapped exactly matches counts of each cub face dims (i.e., each cub face is fully covered exactly)
                            ok_full = True
                            for fcid in range(6):
                                tw,th = face_dims[fcid]
                                cnt = sum(1 for k in occ if k[0]==fcid)
                                if cnt != tw*th:
                                    ok_full = False; break
                            if ok_full:
                                return True
                    # else try next offset
    return False

# ---------- Main solving function ----------
def solve_from_segments(segments):
    W,H,horizontal,vertical = build_walls(segments)
    interior, W, H = find_interior_cells(W,H,horizontal,vertical)
    if not interior:
        print(0); return
    cells_set = set(interior)
    area = len(cells_set)
    triples = find_candidate_triples(area)
    if not triples:
        print(0); return
    # For each triple attempt tiling and folding
    for (A,B,C) in triples:
        # required rectangle areas
        AB = A*B; AC = A*C; BC = B*C
        required_areas = [AB,AB,AC,AC,BC,BC]
        # generate all permutations of assigning which of the 6 faces has which area
        # But also need orientations (w,h) for each assigned area: must be two of A,B,C as sides.
        # We'll generate the multiset of (w,h) for each face by choosing for each area either (x,y) or (y,x) where x,y are corresponding dims
        # For each permutation we will try to tile.
        # Instead of iterating all permutations of areas (720), we will iterate unique permutations to reduce runs.
        from itertools import permutations
        seen_configs = set()
        for perm in set(permutations(required_areas,6)):
            # for each area, get possible (w,h) choices (two possibilities if area corresponds to two different dims, maybe same if square)
            rect_options = []
            possible = True
            for area_val in perm:
                # if area == AB then sides are A,B etc
                if area_val == AB:
                    rect_options.append([(A,B),(B,A)])
                elif area_val == AC:
                    rect_options.append([(A,C),(C,A)])
                elif area_val == BC:
                    rect_options.append([(B,C),(C,B)])
                else:
                    possible = False; break
            if not possible: continue
            # iterate cartesian product of orientation choices (2^6 <= 64)
            for choice in itertools.product(*rect_options):
                key = tuple(choice)
                if key in seen_configs: continue
                seen_configs.add(key)
                # try tile with these rect sizes (order corresponds to perm positions)
                placements = try_tile(cells_set, W, H, list(choice))
                if placements is None:
                    continue
                # Build adjacency
                adj = build_adjacency(placements)
                # Attempt folding simulation
                if try_fold(placements, adj, A,B,C):
                    print(A*B*C)
                    return
    print(0)

# ---------- Entry point ----------
def main():
    data = read_ints()
    if not data:
        return
    it = iter(data)
    N = next(it)
    segments = []
    for _ in range(N):
        x1 = next(it); y1 = next(it); x2 = next(it); y2 = next(it)
        segments.append((x1,y1,x2,y2))
    solve_from_segments(segments)

if __name__ == "__main__":
    main()
