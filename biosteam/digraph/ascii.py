# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2023, Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Plain-text (ASCII) flowsheet diagrams.

Units are drawn inline as ``[ID]`` (the first two characters of their ID in
upper case) and streams as orthogonal lines with arrowheads, for example::

          [B1]<-[H1]<-[F1]
           |           ^
           v     ------|------
          [C1]   |    [A1]   |
           |     |     ^     |
           v     v     |     |
    feed->[S1]->[M1]->[O1]->[R1]->[X1]->[G1]->product
                 ^
                 |
             feed2

Lines are drawn with ``-`` and ``|`` (where lines cross, the vertical line
is drawn on top) and ``+`` is only used for junctions, where streams from
the same source split or streams to the same sink merge (e.g., several feeds
or products of a unit share a manifold). Only feeds and products are
labeled; streams with both a source and a sink are drawn without labels.

The layout is inspired by graphviz's layered drawing:

1. The longest chain of units along the system path (the main line) is
   placed on one row, left to right.
2. Return paths of recycle loops (chains of units that lead back to an
   upstream unit) wrap around the main line, above or below it. Branches
   that feed into or leave the main line are placed on rows above or below.
3. Streams are routed orthogonally with an A* search that penalizes bends
   and crossings, and joins streams with the same source or sink.
4. Whether each loop or branch goes above or below the main line (and the
   spacing) is chosen iteratively by moving to the best neighboring layout
   until a layout repeats. The best layout has the fewest errors (e.g.,
   overlapping lines), then the smallest height.

Unconnected sections of a flowsheet are drawn separately. The layout engine
does not depend on biosteam (see `ascii_flowsheet`).

"""
import heapq

__all__ = (
    'ascii_flowsheet',
    'ascii_from_units',
    'ascii_minimal',
    'ascii_surface',
)

# %% Layout engine (independent of biosteam objects)

RIGHT, DOWN, LEFT, UP = range(4)
DX = (1, 0, -1, 0)
DY = (0, 1, 0, -1)
BEND_COST = 4
CROSS_COST = 1 # Streams may cross (better than going around)
OVERLAP_COST = 500
PORT_BEND_COST = 25 # Bending right next to a port, junction, or arrow (an error)
ADJACENT_COST = 12 # Running right next to another line (an error)
JUNCTION_COST = 2 # Joining lines with the same source or sink
# Extra rows, extra columns, and extra rows between units off the main line
# and the units next to them (they are as close as possible without them)
SPACINGS = tuple([(r, c, d) for r, c in ((0, 0), (1, 0), (2, 0), (3, 0), (0, 1), (1, 1), (2, 1), (2, 2))
                  for d in sorted(set([0, r]))])

FREE, HORIZONTAL, VERTICAL, CROSSING, CORNER, BLOCKED = range(6)
MIXED = object() # Lines from different sources (or to different sinks)


class _Node:
    __slots__ = ('key', 'label', 'order', 'inlets', 'outlets', 'col', 'row',
                 'x', 'y', 'width', 'loop', 'in_side', 'out_side', 'interior', 'free')

    def __init__(self, key, label, order):
        self.key = key
        self.label = label
        self.order = order
        self.inlets = []
        self.outlets = []
        self.loop = False
        self.free = False
        self.in_side = 'L'
        self.out_side = 'R'

    def finalize(self):
        self.interior = max(len(self.label), 1)
        self.width = self.interior + 2

    def text(self):
        return '[' + self.label.center(self.interior) + ']'

    def internal_inlets(self):
        return [i for i in self.inlets if i.src is not None]

    def internal_outlets(self):
        return [i for i in self.outlets if i.dst is not None]


class _Edge:
    __slots__ = ('src', 'src_index', 'label', 'dst_index', 'dst',
                 'src_port', 'dst_port', 'back')

    def __init__(self, src, src_index, label, dst_index, dst):
        self.src = src
        self.src_index = src_index
        self.label = label
        self.dst_index = dst_index
        self.dst = dst
        self.src_port = None
        self.dst_port = None
        self.back = (src is not None and dst is not None and src.order >= dst.order)


def _port_cell(node, port):
    side, i = port
    if side == 'L': return node.x - 1, node.y
    if side == 'R': return node.x + node.width, node.y
    if side == 'T': return node.x + i, node.y - 1
    return node.x + i, node.y + 1

_outward = {'L': LEFT, 'R': RIGHT, 'T': UP, 'B': DOWN}
_inward = {'L': RIGHT, 'R': LEFT, 'T': DOWN, 'B': UP}
_arrow = {'L': '>', 'R': '<', 'T': 'v', 'B': '^'}


class _Grid:
    __slots__ = ('width', 'height', 'chars', 'state', 'collisions',
                 'source', 'sink', 'orientation', 'no_junction', 'cells')

    def __init__(self, width, height):
        self.width = width
        self.height = height
        self.collisions = 0
        self.chars = [[' '] * width for _ in range(height)]
        self.state = [[FREE] * width for _ in range(height)]
        self.source = [[None] * width for _ in range(height)] # Source of lines through cell
        self.sink = [[None] * width for _ in range(height)] # Sink of lines through cell
        self.orientation = [[None] * width for _ in range(height)] # HORIZONTAL, VERTICAL, or CORNER
        self.no_junction = set() # Cells right next to ports
        self.cells = {} # Line cells by source

    def inside(self, x, y):
        return 0 <= x < self.width and 0 <= y < self.height

    def block_text(self, x, y, text):
        if not (0 <= y < self.height):
            self.collisions += len(text)
            return
        chars = self.chars[y]
        state = self.state[y]
        for i, c in enumerate(text):
            xi = x + i
            if 0 <= xi < self.width:
                if state[xi] != FREE and chars[xi] != c: self.collisions += 1
                chars[xi] = c
                state[xi] = BLOCKED
                if c == '|': self.orientation[y][xi] = VERTICAL # Bus of a manifold or stub
            else:
                self.collisions += 1

    def compact(self):
        """Remove columns that only lengthen horizontal lines or widen gaps."""
        chars = self.chars
        W = self.width
        H = self.height
        kept = []
        for x in range(W):
            if kept and x + 1 < W:
                left = kept[-1]
                removable = True
                for y in range(H):
                    row = chars[y]
                    c = row[x]
                    if c == ' ':
                        if row[left] != ' ' and row[x + 1] != ' ': removable = False; break
                    elif c == '-':
                        if row[left] != '-' and row[x + 1] != '-': removable = False; break
                    else:
                        removable = False; break
                if removable: continue
            kept.append(x)
        if len(kept) < W:
            self.chars = chars = [[row[x] for x in kept] for row in chars]
            self.width = W = len(kept)
        # Remove rows that only lengthen vertical lines or widen gaps
        kept = []
        for y in range(H):
            if kept and y + 1 < H:
                above = chars[kept[-1]]
                below = chars[y + 1]
                row = chars[y]
                removable = True
                for x in range(W):
                    c = row[x]
                    if c == ' ':
                        if above[x] != ' ' and below[x] != ' ': removable = False; break
                    elif c == '|':
                        if above[x] != '|' and below[x] != '|': removable = False; break
                    else:
                        removable = False; break
                if removable: continue
            kept.append(y)
        if len(kept) < H:
            self.chars = [chars[y] for y in kept]
            self.height = len(kept)

    def render(self):
        self.compact()
        lines = [''.join(i).rstrip() for i in self.chars]
        while lines and not lines[0]: lines.pop(0)
        while lines and not lines[-1]: lines.pop()
        indent = min([len(i) - len(i.lstrip()) for i in lines if i], default=0)
        return '\n'.join([i[indent:] for i in lines])


def _minimum_bends(x, y, d, gx, gy, gd):
    """Lower bound of the number of bends from (x, y) moving in direction d to
    (gx, gy) arriving in direction gd."""
    along = (gx - x) * DX[d] + (gy - y) * DY[d]
    if d == gd:
        across = (gy - y) if DX[d] else (gx - x)
        if along > 0: return 0 if across == 0 else 2
        return 4
    elif d == (gd + 2) % 4:
        return 2
    else:
        along_goal = (gx - x) * DX[gd] + (gy - y) * DY[gd]
        return 1 if along >= 0 and along_goal > 0 else 3


def _route(grid, starts, goal, goal_dir, goal_open, source, sink):
    """
    A* search for an orthogonal path to the goal port (arriving in the goal
    direction) or to a line going to the same sink (merging at a junction).
    The path may start at any of the given starts (x, y, direction, cost,
    junction, strict); a start with a junction branches off a line from the
    same source. Crossings, bends, and running next to other lines are
    penalized; overlapping other lines is heavily penalized (it is only done
    when there is no other way). Lines are never crossed by lines with the
    same source or sink to avoid confusing crossings with junctions.

    Return the path, the junction it starts from, and the junction it ends at
    (or None).

    """
    state = grid.state
    sources = grid.source
    sinks = grid.sink
    orientation = grid.orientation
    no_junction = grid.no_junction
    W = grid.width
    H = grid.height
    gx, gy = goal
    def h(x, y, d):
        if x == gx and y == gy: return 0
        return abs(x - gx) + abs(y - gy) + BEND_COST * _minimum_bends(x, y, d, gx, gy, goal_dir)
    best = {}
    parents = {}
    initial = {}
    heap = []
    for x, y, d, cost, junction, strict in starts:
        current = (x, y, d)
        if cost < best.get(current, float('inf')):
            best[current] = cost
            parents[current] = None
            initial[current] = (junction, strict)
            heapq.heappush(heap, (cost + h(x, y, d), cost, current))
    while heap:
        f, cost, current = heapq.heappop(heap)
        if best.get(current, None) != cost: continue
        if len(current) == 4: # Merged into a line going to the same sink
            _, mx, my, md = current
            current = parents[current]
            path = []
            while current is not None:
                path.append(current)
                current = parents[current]
            path.reverse()
            return path, initial[path[0]][0], (mx, my, md)
        x, y, d = current
        if x == gx and y == gy:
            if d == goal_dir:
                path = []
                while current is not None:
                    path.append(current)
                    current = parents[current]
                path.reverse()
                return path, initial[path[0]][0], None
            continue
        here = state[y][x]
        for nd in (d, (d + 1) % 4, (d - 1) % 4):
            step = 1
            if nd != d:
                if current in initial: # Leave ports and junctions straight
                    if initial[current][1]: continue
                    step += PORT_BEND_COST
                if here != FREE:
                    if here == BLOCKED: continue
                    step += OVERLAP_COST
                step += BEND_COST
            nx = x + DX[nd]
            ny = y + DY[nd]
            if not (0 <= nx < W and 0 <= ny < H): continue
            horizontal = nd in (LEFT, RIGHT)
            if sink is not None and sinks[ny][nx] == sink:
                o = orientation[ny][nx]
                if ((nx, ny) not in no_junction and
                    (o == CORNER or (o == HORIZONTAL) != horizontal)):
                    new = (None, nx, ny, nd)
                    # Prefer joining close to the arrow
                    new_cost = cost + step + JUNCTION_COST + 2 * (abs(nx - gx) + abs(ny - gy))
                    if nd != d: new_cost += PORT_BEND_COST # Straight line before junctions
                    if new_cost < best.get(new, float('inf')):
                        best[new] = new_cost
                        parents[new] = current
                        heapq.heappush(heap, (new_cost, new_cost, new))
                continue
            if sources[ny][nx] == source: continue
            if nx == gx and ny == gy:
                if not goal_open or nd != goal_dir: continue
                if nd != d: # Straight line before arrows
                    if goal_dir in (LEFT, RIGHT): continue
                    step += PORT_BEND_COST
            else:
                s = state[ny][nx]
                if s == BLOCKED:
                    continue
                elif s != FREE:
                    if (s == VERTICAL and horizontal) or (s == HORIZONTAL and not horizontal):
                        step += CROSS_COST
                    else:
                        step += OVERLAP_COST
            if horizontal:
                for yy in (ny - 1, ny + 1):
                    if 0 <= yy < H and (state[yy][nx] == CROSSING or orientation[yy][nx] in (HORIZONTAL, CORNER)):
                        step += ADJACENT_COST
            else:
                for xx in (nx - 1, nx + 1):
                    if 0 <= xx < W and (state[ny][xx] == CROSSING or orientation[ny][xx] in (VERTICAL, CORNER)):
                        step += ADJACENT_COST
            new = (nx, ny, nd)
            new_cost = cost + step
            if new_cost < best.get(new, float('inf')):
                best[new] = new_cost
                parents[new] = current
                heapq.heappush(heap, (new_cost + h(nx, ny, nd), new_cost, new))
    return None, None, None


def _branches(grid, source):
    """Return starts that branch off lines from the given source."""
    starts = []
    state = grid.state
    orientation = grid.orientation
    no_junction = grid.no_junction
    W = grid.width
    H = grid.height
    for x, y in grid.cells.get(source, ()):
        if (x, y) in no_junction: continue
        o = orientation[y][x]
        for d in range(4):
            if o == HORIZONTAL and d in (LEFT, RIGHT): continue
            if o == VERTICAL and d in (UP, DOWN): continue
            nx = x + DX[d]
            ny = y + DY[d]
            if 0 <= nx < W and 0 <= ny < H and state[ny][nx] == FREE:
                starts.append((nx, ny, d, JUNCTION_COST, (x, y), False))
    return starts


def _draw_path(grid, path, arrow, start_junction, end_junction, source, sink):
    """
    Draw path and return the number of crossings, overlaps, bends, jogs
    (consecutive bends or bends right next to ports, junctions, and arrows),
    and cells right next to parallel lines with a different source and sink
    (lines must have a blank row or column between them).
    """
    chars = grid.chars
    state = grid.state
    sources = grid.source
    sinks = grid.sink
    orientation = grid.orientation
    W = grid.width
    H = grid.height
    cells = grid.cells.setdefault(source, [])
    crossings = overlaps = bends = jogs = adjacent = 0
    own = set([(x, y) for x, y, _ in path])
    def next_to_line(x, y, horizontal):
        if horizontal:
            neighbors = ((x, y - 1), (x, y + 1))
            kinds = (HORIZONTAL, CORNER, CROSSING)
        else:
            neighbors = ((x - 1, y), (x + 1, y))
            kinds = (VERTICAL, CORNER, CROSSING)
        for nx, ny in neighbors:
            if (0 <= nx < W and 0 <= ny < H and (nx, ny) not in own
                and (state[ny][nx] == CROSSING or orientation[ny][nx] in kinds)
                and sources[ny][nx] != source and sinks[ny][nx] != sink):
                return True
        return False
    def claim(x, y):
        sources[y][x] = source if sources[y][x] in (None, source) else MIXED
        sinks[y][x] = sink if sinks[y][x] in (None, sink) else MIXED
        cells.append((x, y))
    if start_junction:
        x, y = start_junction
        chars[y][x] = '+'
        orientation[y][x] = CORNER
        if state[y][x] != BLOCKED: state[y][x] = CORNER
        sinks[y][x] = MIXED # Lines go to different sinks
    if end_junction:
        x, y, d = end_junction
        chars[y][x] = '+'
        orientation[y][x] = CORNER
        if state[y][x] != BLOCKED: state[y][x] = CORNER
        sources[y][x] = MIXED # Lines come from different sources
        path = path + [end_junction]
    n = len(path)
    for i, (x, y, d) in enumerate(path):
        if end_junction and i == n - 1: break
        s = state[y][x]
        if i == n - 1:
            chars[y][x] = arrow
            state[y][x] = BLOCKED
            continue
        nd = path[i + 1][2]
        if nd != d:
            if s != FREE: overlaps += 1
            if i == 0 or i == n - 2 or path[i - 1][2] != d: jogs += 1
            bends += 1
            chars[y][x] = '-' # '+' is reserved for junctions
            state[y][x] = orientation[y][x] = CORNER
        elif d in (LEFT, RIGHT):
            if s == VERTICAL:
                crossings += 1
                state[y][x] = CROSSING
            else:
                if s != FREE: overlaps += 1
                chars[y][x] = '-'
                state[y][x] = orientation[y][x] = HORIZONTAL
        else:
            if s == HORIZONTAL:
                crossings += 1
                state[y][x] = CROSSING
            else:
                if s != FREE: overlaps += 1
                state[y][x] = orientation[y][x] = VERTICAL
            chars[y][x] = '|'
        if state[y][x] != CROSSING:
            if next_to_line(x, y, d in (LEFT, RIGHT)): adjacent += 1
            claim(x, y)
    return crossings, overlaps, bends, jogs, adjacent


class _Layout:

    def __init__(self, nodes, edges):
        self.nodes = nodes
        self.edges = edges
        self.internal = [e for e in edges if e.src is not None and e.dst is not None]
        self.find_loops()
        self.assign_columns()
        self.find_main_line()
        self.find_groups()

    def find_main_line(self):
        if not self.internal:
            self.main = []
            return
        nodes = [n for n in self.nodes if n not in self.loop_nodes]
        best = {}
        parent = {}
        for n in nodes:
            best[n] = 1
            parent[n] = None
            for e in n.internal_inlets():
                if not e.back and e.src in best and best[e.src] + 1 > best[n]:
                    best[n] = best[e.src] + 1
                    parent[n] = e.src
        end = max(nodes, key=lambda n: (best[n], -n.order))
        main = []
        while end is not None:
            main.append(end)
            end = parent[end]
        self.main = main[::-1]

    def find_loops(self):
        """Find return paths of recycle loops (chains of units with one
        internal inlet and one internal outlet that end in a recycle)."""
        used = set()
        loops = []
        for e in sorted(self.internal, key=lambda e: e.src.order):
            if not e.back: continue
            target = e.dst
            chain = []
            current = e.src
            while True:
                if current in used or current is target: break
                ins = current.internal_inlets()
                outs = current.internal_outlets()
                if len(ins) != 1 or len(outs) != 1: break
                if len(outs) != len(current.outlets): break # Products stay on the main line
                if outs[0] is not e and outs[0].dst is not (chain[0] if chain else None): break
                chain.insert(0, current)
                previous = ins[0].src
                if previous in chain or ins[0].back: break
                current = previous
            if not chain: continue
            origin = chain[0].internal_inlets()[0].src
            if origin in chain or origin is target: continue
            loops.append((origin, chain, target))
            used.update(chain)
        self.loops = loops

    def assign_columns(self):
        nodes = self.nodes
        loop_nodes = set()
        for _, chain, _ in self.loops: loop_nodes.update(chain)
        for _ in range(2):
            others = [n for n in nodes if n not in loop_nodes]
            for n in others:
                n.col = max([e.src.col + 1 for e in n.internal_inlets()
                             if not e.back and e.src not in loop_nodes], default=0)
            for n in reversed(others): # Pull nodes toward their successors
                cols = [e.dst.col for e in n.internal_outlets()
                        if not e.back and e.dst not in loop_nodes]
                if cols:
                    col = min(cols) - 1
                    if col > n.col: n.col = col
            # Loop return paths wrap around the main line from the origin
            # back to the target (placed in `place_rows`)
            loops = []
            for origin, chain, target in self.loops:
                if origin.col > target.col:
                    for n in chain: n.loop = True
                    loops.append((origin, chain, target))
                else:
                    for n in chain:
                        n.loop = False
                        loop_nodes.discard(n)
            if len(loops) == len(self.loops): break
            self.loops = loops
        self.loop_nodes = loop_nodes

    def find_groups(self):
        """Groups of units placed together above or below the main line."""
        main = set(self.main)
        groups = [list(chain) for _, chain, _ in self.loops]
        rest = [n for n in self.nodes if n not in main and n not in self.loop_nodes]
        free = [n for n in rest if not (n.internal_inlets() or n.internal_outlets())]
        for n in free: n.free = True
        self.free_nodes = free
        rest = [n for n in rest if not n.free]
        rest_set = set(rest)
        seen = set()
        for n in rest:
            if n in seen: continue
            group = []
            stack = [n]
            seen.add(n)
            while stack:
                m = stack.pop()
                group.append(m)
                for e in m.internal_inlets() + m.internal_outlets():
                    other = e.src if e.dst is m else e.dst
                    if other in rest_set and other not in seen:
                        seen.add(other)
                        stack.append(other)
            group.sort(key=lambda n: n.order)
            groups.append(group)
        self.groups = groups
        self.n_loops = len(self.loops)

    def place_rows(self, signs):
        occupied = {}
        for n in self.nodes: n.row = None
        for n in self.main:
            n.row = 0
            occupied[n.col, 0] = n
        def free(col, row, width=0):
            return all([(c, row) not in occupied for c in range(col - width, col + width + 1)])
        n_loops = self.n_loops
        for index, (group, sign) in enumerate(zip(self.groups, signs)):
            if not group or index < n_loops: continue
            connected = [n for n in group if n.internal_inlets() or n.internal_outlets()]
            unconnected = [n for n in group if n not in connected]
            placed = set()
            pending = list(connected)
            for direction in (1, -1):
                remaining = []
                for n in (pending if direction == 1 else reversed(pending)):
                    if direction == 1:
                        rows = [e.src.row for e in n.internal_inlets() if e.src.row is not None]
                    else:
                        rows = [e.dst.row for e in n.internal_outlets() if e.dst.row is not None]
                    if not rows and direction == 1:
                        remaining.append(n)
                        continue
                    row = rows[0] if rows else 0
                    if row == 0 or (row > 0) != (sign > 0):
                        row = 0 if free(n.col, 0) else sign
                    while not free(n.col, row): row += sign
                    n.row = row
                    occupied[n.col, row] = n
                    placed.add(n)
                pending = remaining
            if unconnected:
                rows = [r for (_, r) in occupied]
                row = (max(rows) if sign > 0 else min(rows)) + sign
                for col, n in enumerate(unconnected):
                    n.col = col
                    n.row = row
                    occupied[col, row] = n
        for loop, sign in zip(self.loops, signs):
            self.place_loop(*loop, sign, occupied, free)
        if self.free_nodes: # Units without connections to other units
            row = max([r for (_, r) in occupied]) + 1 if occupied else 0
            for n in self.free_nodes: n.row = row

    def place_loop(self, origin, chain, target, sign, occupied, free):
        """
        Place the return path of a recycle loop on a U-shaped sequence of
        cells: away from the main line at the origin's column, back along
        a row toward the target's column, and into the target. Units are
        placed as close to the origin (the source of the loop) as possible.
        """
        X = origin.col
        Y = target.col
        m = len(chain)
        h = max(1, -(-(m - (X - Y - 1)) // 2))
        if sign < 0:
            r0 = min(origin.row, target.row)
        else:
            r0 = max(origin.row, target.row)
        rows = [r for (_, r) in occupied]
        limit = max(rows) - min(rows) + h + 3
        cells = None
        # Find the closest free cells, preferably with a clear line from the
        # origin and to the target (rows beyond the drawing are always free)
        for clear_approach in (True, False):
            base = r0 + sign
            for _ in range(limit):
                top = base + sign * (h - 1)
                candidate = (
                    [(X, base + sign * k) for k in range(h)]
                    + [(c, top) for c in range(X - 1, Y, -1)]
                    + [(Y, base + sign * k) for k in reversed(range(h))]
                )
                approach = [(X, r) for r in range(origin.row + sign, base, sign)]
                approach += [(Y, r) for r in range(target.row + sign, base, sign)]
                if (all([free(c, r) for c, r in candidate])
                    and (not clear_approach or all([free(c, r) for c, r in approach]))):
                    cells = candidate
                    break
                base += sign
            if cells: break
        cells = cells[:m] # Units are placed as close to the source as possible
        for n, (c, r) in zip(chain, cells):
            n.col = c
            n.row = r
            occupied[c, r] = n
        def facing(n, other):
            dx = other.col - n.col
            dy = other.row - n.row
            if dx == 0: return 'B' if dy > 0 else 'T'
            return 'R' if dx > 0 else 'L'
        neighbors = [origin, *chain, target]
        for i, n in enumerate(chain, 1):
            n.in_side = facing(n, neighbors[i - 1])
            n.out_side = facing(n, neighbors[i + 1])

    def assign_ports(self):
        for n in self.nodes:
            vertical = 'T' if n.row < 0 else 'B'
            remaining = []
            if n.loop: # Sides follow the loop's return path
                for edges, side, attr in ((n.inlets, n.in_side, 'dst_port'),
                                          (n.outlets, n.out_side, 'src_port')):
                    for e in edges:
                        internal = e.src is not None and e.dst is not None
                        if internal and side in ('L', 'R'):
                            setattr(e, attr, (side, 0))
                        else:
                            remaining.append((e, attr))
            else:
                vertical_inlet = self.vertical_inlet(n)
                def in_score(e):
                    if e.src is None: return 1
                    if e.back or e.src.loop or e is vertical_inlet: return None
                    if e.src.row == n.row and e.src.col < n.col: return 0
                    return 2
                def out_score(e):
                    if e.dst is None: return 1
                    if e.back or e.dst.loop: return None
                    if e.dst.row == n.row and e.dst.col > n.col: return 0
                    return 2
                forward = [e for e in n.outlets if out_score(e) is not None]
                if len(forward) > 1: # All forward outlets leave through a manifold
                    for e in forward: e.src_port = ('R', 0)
                    # Recycles also leave through the manifold (from an end of the bus)
                    recycles = [e for e in n.outlets if e not in forward
                                and e.dst is not None and e.back and not e.dst.loop]
                    for e in recycles: e.src_port = ('R', 0)
                    remaining.extend([(e, 'src_port') for e in n.outlets
                                      if e not in forward and e not in recycles])
                    sides = ((n.inlets, in_score, 'L', 'dst_port'),)
                else:
                    sides = ((n.inlets, in_score, 'L', 'dst_port'),
                             (n.outlets, out_score, 'R', 'src_port'))
                for edges, score, side, attr in sides:
                    scored = [(score(e), i, e) for i, e in enumerate(edges)]
                    scored = [i for i in scored if i[0] is not None]
                    if scored:
                        _, _, best = min(scored, key=lambda x: (x[0], x[1]))
                        if best.src is None or best.dst is None:
                            # All feeds (or products) share the side through a manifold
                            chosen = [e for e in edges if e.src is None or e.dst is None]
                        elif attr == 'dst_port':
                            # Inlets from upstream merge before the arrow
                            chosen = [e for e in edges if score(e) in (0, 2)]
                        else:
                            chosen = [best]
                        for e in chosen: setattr(e, attr, (side, 0))
                    else:
                        chosen = []
                    for e in edges:
                        if e not in chosen: remaining.append((e, attr))
            top = []
            bottom = []
            for e, attr in remaining:
                other = e.src if attr == 'dst_port' else e.dst
                if n.loop and other is not None:
                    side = n.in_side if attr == 'dst_port' else n.out_side
                elif other is None:
                    side = vertical
                elif e.back and attr == 'src_port':
                    # Recycles leave away from the main line
                    side = 'T' if n.row <= 0 else 'B'
                elif other.row < n.row:
                    side = 'T'
                elif other.row > n.row:
                    side = 'B'
                elif e.back:
                    side = 'T'
                else:
                    side = vertical
                (top if side == 'T' else bottom).append((e, attr))
            # Feeds (and products) on the same side share a single stub
            w = n.width
            slots = {}
            for side, items in (('T', top), ('B', bottom)):
                feeds = [i for i in items if i[0].src is None]
                products = [i for i in items if i[0].dst is None]
                internal = [i for i in items if i[0].src is not None and i[0].dst is not None]
                def key(item):
                    e, attr = item
                    other = e.src if attr == 'dst_port' else e.dst
                    return (other.col, other.row)
                internal.sort(key=key)
                inlets = [i for i in internal if i[1] == 'dst_port']
                outlets = [i for i in internal if i[1] == 'src_port']
                # Lines to (or from) the same side share an anchor and are joined at junctions
                slots[side] = [i for i in (feeds, inlets, outlets, products) if i]
            for side, other in (('T', 'B'), ('B', 'T')):
                while len(slots[side]) > w and len(slots[other]) < w:
                    internal = [i for i in slots[side] if i[0][0].src is not None and i[0][0].dst is not None]
                    if not internal: break
                    slot = internal[-1]
                    slots[side].remove(slot)
                    slots[other].insert(-1 if slots[other] and slots[other][-1][0][0].dst is None else len(slots[other]), slot)
            w = max(w, len(slots['T']), len(slots['B'])) # Only widen in extreme cases
            for side, items in slots.items():
                k = len(items)
                if k == 0: continue
                stubs = any([e.src is None or e.dst is None for slot in items for e, _ in slot])
                if k == 1:
                    offsets = [1]
                elif k == 2:
                    offsets = [1, w - 1]
                else:
                    offsets = [round(j * (w - 1) / (k - 1)) for j in range(k)]
                for slot, offset in zip(items, offsets):
                    for e, attr in slot: setattr(e, attr, (side, offset))
            n.interior = w - 2
            n.width = w

    def stubs(self, node):
        """Return feeds/products by port."""
        stubs = {}
        for e in node.inlets + node.outlets:
            if e.src is not None and e.dst is not None: continue
            port = e.dst_port if e.src is None else e.src_port
            stubs.setdefault(port, []).append(e)
        return stubs

    def has_bus(self, node):
        """Return whether the outlets of a unit leave through a manifold."""
        if node.loop: return False
        forward = [e for e in node.outlets
                   if e.dst is None or not (e.back or e.dst.loop)]
        return len(forward) > 1

    def vertical_inlet(self, node):
        """
        Return the inlet of a unit that comes straight up from the manifold bus
        of a unit below it (into the bottom of the unit), if any.
        """
        if node.loop or node.free: return None
        inlets = node.internal_inlets()
        if len(inlets) != 1: return None
        e = inlets[0]
        src = e.src
        if e.back or src.loop or node.row >= src.row or not self.has_bus(src): return None
        lo, hi = sorted([src.row, node.row])
        if any([n is not node and n is not src and not n.free and n.col == src.col
                and lo <= n.row <= hi for n in self.nodes]): return None
        return e

    def manifold(self, node):
        """
        Return the branches of the right manifold of a node as a tuple of
        edges aligned with the row of their target, offsets (relative to
        the node row) of other branches (products and edges to the top or
        bottom of their target), and recycles (which leave from an end of the
        bus). Return None if the node has no manifold.
        """
        edges = [e for e in node.outlets if e.src_port == ('R', 0)]
        recycles = [e for e in edges if e.dst is not None and e.back]
        edges = [e for e in edges if e not in recycles]
        if len(edges) < 2: return None
        bus = node.x + node.width if node.x is not None else None
        vertical = [e for e in edges if e.dst is not None and e.dst_port[0] in ('T', 'B')
                    and e is self.vertical_inlet(e.dst)
                    and (e.dst.x is None or bus is None or e.dst.x + e.dst_port[1] == bus)]
        edges = [e for e in edges if e not in vertical]
        aligned = []
        rows = set()
        for e in sorted([e for e in edges if e.dst is not None and e.dst_port[0] == 'L'],
                        key=lambda e: e.dst.col):
            if e.dst.row in rows: continue # Only one branch per row can be aligned
            if e.dst.row != node.row:
                # The bus cannot pass other units in the same column
                lo, hi = sorted([node.row, e.dst.row])
                if any([n is not node and not n.free and n.col == node.col and lo <= n.row <= hi
                        for n in self.nodes]): continue
            rows.add(e.dst.row)
            aligned.append(e)
        others = [e for e in edges if e not in aligned]
        used = set([0]) if any([e.dst.row == node.row for e in aligned]) else set()
        away = -1 if node.row < 0 else 1 # Away from the main line
        offsets = {}
        # Branches to units go toward their unit; products go away from the main line
        others.sort(key=lambda e: e.dst is None)
        for e in others:
            if e.dst is None:
                sign = away
            elif e.dst.row == node.row: # Opposite to products
                sign = -away
            else:
                sign = 1 if e.dst.row > node.row else -1
            offset = 0
            while offset in used: offset += 2 * sign
            offsets[e] = offset
            used.add(offset)
        return aligned, offsets, recycles, vertical

    def row_extent(self, node):
        """Horizontal extent of a node in its row, including feeds and products
        on its sides and its top band."""
        left = right = 1
        for (side, _), edges in self.stubs(node).items():
            length = max([len(e.label) for e in edges]) + 3
            if side == 'L': left = max(left, length + 1)
            elif side == 'R': right = max(right, length + 1)
        x0 = node.x - left
        x1 = node.x + node.width + right
        extents = self.band_extents(node)
        if 'T' in extents:
            lo, hi = extents['T']
            x0 = min(x0, node.x + lo)
            x1 = max(x1, node.x + hi)
        return x0, x1

    def uses_side(self, node, side):
        """Return whether streams between units connect to the given side."""
        return (
            any([e.src is not None and e.dst_port[0] == side for e in node.inlets])
            or any([e.dst is not None and e.src_port[0] == side for e in node.outlets])
        )

    def band_rows(self, node, side):
        """Number of rows taken above (T) or below (B) the node by feeds and
        products."""
        rows = 0
        for (s, _), edges in self.stubs(node).items():
            N = len(edges)
            if s == side:
                rows = max(rows, 2 * N + 1)
            elif s == 'L' and N > 1: # Manifold (feeds every other row, starting at the node row)
                above = 2 * ((N - 1) // 2)
                rows = max(rows, above if side == 'T' else 2 * (N - 1) - above)
        manifold = self.manifold(node)
        if manifold:
            sign = -1 if side == 'T' else 1
            for offset in manifold[1].values():
                if offset * sign > 0: rows = max(rows, abs(offset))
        return rows

    def band_extents(self, node):
        """Horizontal extent of feeds and products drawn above (T) or below (B)
        the node, relative to the node."""
        extents = {}
        def add(side, lo, hi):
            if side in extents:
                lo0, hi0 = extents[side]
                lo = min(lo, lo0)
                hi = max(hi, hi0)
            extents[side] = (lo, hi)
        for (side, x), edges in self.stubs(node).items():
            N = len(edges)
            length = max([len(e.label) for e in edges])
            feed = edges[0].src is None
            if side in ('T', 'B'):
                if N == 1:
                    lo, hi = (x - length + 1, x) if feed else (x, x + length - 1)
                else:
                    lo, hi = (x - length - 1, x) if feed else (x, x + length + 2)
                add(side, lo, hi)
            elif side == 'L' and N > 1: # Manifold
                add('T', -length - 3, -2)
                add('B', -length - 3, -2)
        manifold = self.manifold(node)
        if manifold:
            w = node.width
            for e, offset in manifold[1].items():
                if offset == 0: continue
                side = 'T' if offset < 0 else 'B'
                add(side, w, w + 2 + len(e.label) if e.dst is None else w + 1)
            for e in manifold[0]: # The bus passes between rows to aligned branches
                if e.dst.row != node.row:
                    add('B' if e.dst.row > node.row else 'T', w, w)
        return extents

    def position(self, extra):
        extra_rows, extra_cols, extra_decoupled = extra
        nodes = self.nodes
        placed = [n for n in nodes if not n.free]
        cols = sorted(set([n.col for n in placed]))
        rows = sorted(set([n.row for n in nodes]))
        # Margin for recycles going around everything (two rows each and room for arrows)
        n_back = 2 * len([e for e in self.internal if e.back]) + 3
        # Horizontal positions: units are aligned by column, but are shifted
        # right when needed by labels in their row or by connected units
        row_end = {r: n_back for r in rows} # End of last unit or label in row
        row_last = {r: None for r in rows} # Edge leaving right side of last unit in row
        stub_end = {} # End of last top/bottom stub label by row and side
        def side_edge(n, side):
            for e in n.inlets + n.outlets:
                if (e.dst_port if e.dst is n else e.src_port) == (side, 0): return e
        def place(n, x):
            left = side_edge(n, 'L')
            last = row_last[n.row]
            if left is not None and (left.src is None or left.dst is None):
                feeds = [e for e in n.inlets if e.dst_port == ('L', 0)]
                gap = max([len(e.label) for e in feeds]) + (3 if len(feeds) == 1 else 4)
                if last is not None: gap = max(gap, 5)
            elif left is not None and left is last:
                # Room for junctions before the arrow
                gap = 3 if len([e for e in n.inlets if e.dst_port == ('L', 0)]) > 1 else 2
            elif left is not None or last is not None:
                gap = 5
            else:
                gap = 2
            x = max(x, row_end[n.row] + gap)
            extents = self.band_extents(n)
            for side, (lo, hi) in extents.items():
                x = max(x, stub_end.get((n.row, side), -2) + 3 - lo)
            n.x = x
            right = side_edge(n, 'R')
            end = x + n.width
            manifold = self.manifold(n)
            if manifold:
                aligned, offsets, _, _ = manifold
                end += 1
                right = None
                for e in aligned:
                    if e.dst.row == n.row: right = e
                for e, offset in offsets.items():
                    if offset == 0: end += 2 + (len(e.label) if e.dst is None else 0)
            elif right is not None and (right.src is None or right.dst is None):
                products = [e for e in n.outlets if e.src_port == ('R', 0)]
                end += max([len(e.label) for e in products]) + (2 if len(products) == 1 else 3)
                right = None
            row_end[n.row] = end
            row_last[n.row] = right
            for side, (lo, hi) in extents.items():
                stub_end[n.row, side] = x + hi
        for n in nodes: n.x = None
        base = n_back
        previous = None
        for c in cols:
            if previous is not None:
                base = min([n.x for n in placed if n.col == previous]) + 1 + extra_cols # Keep columns in order
            members = [n for n in placed if n.col == c]
            members.sort(key=lambda n: (n.loop, abs(n.row)))
            for n in members:
                x = base
                vertical = self.vertical_inlet(n)
                for e in n.internal_inlets():
                    p = e.src
                    if e is vertical: continue
                    if not e.back and p.x is not None and p.col < c:
                        x = max(x, p.x + p.width + (3 if e.src_port == ('R', 0) and self.manifold(p) else 2))
                if n.loop: # Align with units directly above or below in the loop
                    for side, edges, attr in ((n.in_side, n.internal_inlets(), 'src'),
                                              (n.out_side, n.internal_outlets(), 'dst')):
                        if side in ('T', 'B'):
                            for e in edges:
                                other = getattr(e, attr)
                                if other.x is not None and other.col == c: x = max(x, other.x)
                vertical = self.vertical_inlet(n)
                if vertical and vertical.src.x is not None:
                    x = max(x, vertical.src.x + vertical.src.width - vertical.dst_port[1])
                place(n, x)
            previous = c
        # Pull units off the main line toward the units they feed to keep lines short
        main = set(self.main)
        for n in sorted(placed, key=lambda n: -n.x):
            if n in main or n.loop or self.manifold(n) or self.vertical_inlet(n): continue
            outlets = n.internal_outlets()
            if not outlets or len(outlets) != len(n.outlets): continue
            if any([e.back or e.dst.x is None or e.dst.x <= n.x for e in outlets]): continue
            desired = []
            for e in outlets:
                side, offset = e.dst_port
                if side == 'L':
                    merged = len([i for i in e.dst.inlets if i.dst_port == ('L', 0)]) > 1
                    gap = 2 if e.dst.row == n.row and not merged else 3
                    desired.append(e.dst.x - n.width - gap)
                else:
                    desired.append(e.dst.x + offset - n.width - 2)
            x = min(desired)
            successors = set([e.dst for e in outlets])
            for m in placed: # Do not pass other units in the same row
                if m is not n and m not in successors and m.row == n.row and m.x > n.x:
                    x = min(x, m.x - n.width - 5)
            if x > n.x:
                if self.band_extents(n): continue # Feeds or products above or below
                n.x = x
        for n in self.free_nodes: place(n, n_back)
        width = max(list(row_end.values()) + list(stub_end.values())) + n_back
        # Vertical positions
        row_y = {}
        top = max([self.band_rows(n, 'T') for n in nodes if n.row == rows[0]])
        y = n_back + top
        previous = None
        for r in rows:
            if previous is not None:
                upper = [n for n in nodes if n.row == previous]
                lower = [n for n in nodes if n.row == r]
                # Lines to the top or bottom of units need room for the arrow
                vertical = (
                    any([self.uses_side(n, 'B') for n in upper])
                    or any([self.uses_side(n, 'T') for n in lower])
                )
                gap = 3 if vertical else 2
                # Feeds and products below a unit only push down units they overlap with
                for b in lower:
                    above = self.band_rows(b, 'T')
                    gap = max(gap, above + 1)
                    b0, b1 = self.row_extent(b)
                    for a in upper:
                        below = self.band_rows(a, 'B')
                        if not below: continue
                        lo, hi = self.band_extents(a)['B']
                        if a.x + lo <= b1 + 1 and b0 - 1 <= a.x + hi:
                            gap = max(gap, below + above + 1)
                gap += getattr(self, 'row_extra', {}).get(r, extra_rows)
                y += gap
            row_y[r] = y
            previous = r
        height = y + 1 + max([self.band_rows(n, 'B') for n in nodes if n.row == rows[-1]]) + n_back
        for n in nodes: n.y = row_y[n.row]
        self.decouple_rows(extra_decoupled)
        return width, height

    def decouple_rows(self, extra_rows):
        """
        Move units off the main line as close to it as the units between them
        allow (units in a row do not need to be at the same height). Units
        connected side by side move together.
        """
        nodes = [n for n in self.nodes if not n.free]
        extents = {n: self.row_extent(n) for n in nodes}
        def overlap(a, b):
            a0, a1 = extents[a]
            b0, b1 = extents[b]
            return a0 <= b1 + 1 and b0 <= a1 + 1
        tight = getattr(self, 'tight', set())
        self.decoupled = decoupled = []
        def gap(upper, lower, extra):
            vertical = self.uses_side(upper, 'B') or self.uses_side(lower, 'T')
            bands = self.band_rows(upper, 'B') + self.band_rows(lower, 'T') + 1
            return max(3 if vertical else 2, bands) + extra
        for sign in (-1, 1):
            rows = sorted(set([n.row for n in nodes if n.row * sign > 0]), key=abs)
            for r in rows:
                members = [n for n in nodes if n.row == r]
                # Units connected side by side in the row move together
                groups = []
                seen = set()
                for n in members:
                    if n in seen: continue
                    group = []
                    stack = [n]
                    seen.add(n)
                    while stack:
                        m = stack.pop()
                        group.append(m)
                        for e in m.internal_inlets() + m.internal_outlets():
                            other = e.src if e.dst is m else e.dst
                            port = e.dst_port if e.dst is m else e.src_port
                            if (other in members and other not in seen
                                and port[0] in ('L', 'R') and other.row == r):
                                seen.add(other)
                                stack.append(other)
                    groups.append(group)
                for group in groups:
                    key = min([n.order for n in group])
                    decoupled.append(key)
                    extra = 0 if key in tight else extra_rows
                    ys = []
                    for n in group:
                        for m in nodes:
                            if m in group or not (0 <= m.row * sign < abs(r)): continue
                            if not overlap(n, m): continue
                            if sign < 0:
                                ys.append(m.y - gap(n, m, extra))
                            else:
                                ys.append(m.y + gap(m, n, extra))
                    if not ys: continue
                    y = min(ys) if sign < 0 else max(ys)
                    if (y - group[0].y) * sign < 0: # Only move toward the main line
                        for n in group: n.y = y

    def draw(self, extra):
        width, height = self.position(extra)
        grid = _Grid(width, height)
        nodes = self.nodes
        for n in nodes: grid.block_text(n.x, n.y, n.text())
        labels_drawn = {True: [], False: []} # Rows and extents of feeds (True) and products (False)
        def draw_label(feed, x, y, text):
            grid.block_text(x, y, text)
            labels_drawn[feed].append((y, x, x + len(text) - 1))
            # Keep a space between labels and lines
            for space in ((x - 1, x - 2) if feed else (x + len(text), x + len(text) + 1)):
                if grid.inside(space, y) and grid.state[y][space] == FREE:
                    grid.state[y][space] = BLOCKED
        # Feeds and products
        for n in nodes:
            for (side, i), edges in self.stubs(n).items():
                feed = edges[0].src is None
                labels = [e.label for e in edges]
                N = len(labels)
                x, y = _port_cell(n, (side, i))
                if side == 'R' and self.manifold(n): continue
                if side in ('L', 'R'):
                    if N == 1:
                        if side == 'L':
                            text = labels[0] + '->'
                            draw_label(True, x - len(text) + 1, y, text)
                        else:
                            draw_label(False, x, y, '->' + labels[0])
                        continue
                    # Manifold
                    grid.block_text(x, y, '>')
                    x -= 1
                    rows = [y + 2 * k - 2 * ((N - 1) // 2) for k in range(N)]
                    for yy in range(rows[0], rows[-1] + 1):
                        if yy != y and yy not in rows: grid.block_text(x, yy, '|')
                    if y not in rows: grid.block_text(x, y, '+')
                    for label, yy in zip(labels, rows):
                        text = label + '-+'
                        draw_label(True, x - len(text) + 1, yy, text)
                    continue
                sign = -1 if side == 'T' else 1
                if feed:
                    grid.block_text(x, y, _arrow[side])
                    grid.block_text(x, y + sign, '|')
                else:
                    grid.block_text(x, y, '|')
                    grid.block_text(x, y + sign, '|' if N > 1 else ('v' if side == 'B' else '^'))
                if N == 1:
                    label = labels[0]
                    draw_label(feed, x - len(label) + 1 if feed else x, y + 2 * sign, label)
                    continue
                # Comb
                for k, label in enumerate(labels):
                    yy = y + sign * (2 + 2 * k)
                    if k: grid.block_text(x, yy - sign, '|')
                    if feed:
                        text = label + '-+'
                        draw_label(True, x - len(text) + 1, yy, text)
                    else:
                        draw_label(False, x, yy, '+->' + label)
        # Right manifolds
        drawn = set() # Streams drawn as part of a manifold
        starts = {} # Start cell, direction, and source of lines from manifolds
        for n in nodes:
            manifold = self.manifold(n)
            if not manifold: continue
            aligned, offsets, recycles, vertical = manifold
            offsets = offsets.copy()
            for e in aligned: offsets[e] = _port_cell(e.dst, e.dst_port)[1] - n.y
            x = n.x + n.width
            ends = [_port_cell(e.dst, e.dst_port)[1] - n.y for e in vertical]
            lo = min([0] + list(offsets.values()) + ends)
            hi = max([0] + list(offsets.values()) + ends)
            junctions = set([n.y + i for i in offsets.values()] + [n.y + i for i in ends])
            for e, end in zip(vertical, ends): # The bus goes straight into the unit
                grid.block_text(x, n.y + end, _arrow[e.dst_port[0]])
                drawn.add(e)
            for yy in range(n.y + lo, n.y + hi + 1):
                if yy != n.y and yy not in junctions: grid.block_text(x, yy, '|')
            if n.y not in junctions: grid.block_text(x, n.y, '+')
            up = any([not (e.dst_port[0] == 'B' or e.dst.row > n.row) for e in recycles])
            down = any([e.dst_port[0] == 'B' or e.dst.row > n.row for e in recycles])
            for e, offset in offsets.items():
                if e.dst is None:
                    draw_label(False, x, n.y + offset, '+->' + e.label)
                else:
                    # A line to a unit leaving an end of the bus is a corner, not a junction
                    corner = (offset == lo < 0 and not up) or (offset == hi > 0 and not down)
                    grid.block_text(x, n.y + offset, '-' if corner else '+')
                    cell = (x + 1, n.y + offset)
                    starts[e] = (cell, RIGHT, (n, cell))
            for e in recycles: # Continue the bus up (or down) to the recycle
                if e.dst_port[0] == 'B' or e.dst.row > n.row:
                    starts[e] = ((x, n.y + hi + 1), DOWN, (n, DOWN))
                else:
                    starts[e] = ((x, n.y + lo - 1), UP, (n, UP))
        # Reserve ports
        ports = []
        for e in self.internal:
            if e in drawn: continue
            start = starts[e][0] if e in starts else _port_cell(e.src, e.src_port)
            goal = _port_cell(e.dst, e.dst_port)
            ports.append((e, start, goal))
            for x, y in (start, goal):
                if grid.inside(x, y) and grid.state[y][x] == FREE: grid.state[y][x] = BLOCKED
            grid.no_junction.add(start)
        # Shorter lines first; recycles to units further left go outside
        ports.sort(key=lambda x: (abs(x[1][0] - x[2][0]) + abs(x[1][1] - x[2][1]),
                                  x[0].back, -x[2][0] if x[0].back else 0))
        crossings = overlaps = bends = jogs = adjacent = failed = 0
        paths = {}
        used = set()
        connected = set() # Sources and sinks already connected by a line
        for e, start, goal in ports:
            if not (grid.inside(*start) and grid.inside(*goal)):
                failed += 1
                continue
            if e in starts:
                _, start_dir, source = starts[e]
            else:
                start_dir = _outward[e.src_port[0]]
                source = (e.src, e.src_port)
            sink = (e.dst, e.dst_port)
            if (source, sink) in connected: continue # Same line as another stream
            connected.add((source, sink))
            goal_dir = _inward[e.dst_port[0]]
            options = _branches(grid, source)
            if start not in used:
                options.append((*start, start_dir, 0, None, start_dir in (LEFT, RIGHT)))
                grid.state[start[1]][start[0]] = FREE
            goal_open = goal not in used
            if goal_open: grid.state[goal[1]][goal[0]] = FREE
            path, start_junction, end_junction = _route(
                grid, options, goal, goal_dir, goal_open, source, sink
            )
            if path is not None:
                if start_junction is None: used.add(start)
                if end_junction is None: used.add(goal)
            for x, y in (start, goal): # Keep reserving ports that are not used yet
                if (x, y) not in used and grid.state[y][x] == FREE: grid.state[y][x] = BLOCKED
            if path is None:
                failed += 1
                continue
            c, o, b, j, a = _draw_path(grid, path, _arrow[e.dst_port[0]],
                                       start_junction, end_junction, source, sink)
            crossings += c
            overlaps += o
            bends += b
            jogs += j
            adjacent += a
            paths[e] = [(x, y) for x, y, _ in path]
            # Keep other lines from crossing right next to ports and arrows
            if start_junction is None and e.src_port[0] in ('L', 'R'):
                x, y, _ = path[0]
                grid.state[y][x] = BLOCKED
            if end_junction is None and e.dst_port[0] in ('L', 'R') and len(path) > 1:
                x, y, _ = path[-2]
                grid.state[y][x] = BLOCKED
        outside = 0
        for origin, chain, target in self.loops:
            # Lines of other streams stay inside the loop unless they enclose it
            above = chain[0].row < origin.row
            outer = min([n.y for n in chain]) if above else max([n.y for n in chain])
            lo = min([target.x] + [n.x for n in chain])
            hi = max([origin.x + origin.width] + [n.x + n.width for n in chain]) - 1
            members = set(chain)
            for e, cells in paths.items():
                if e.src in members or e.dst in members: continue
                xs = [x for x, _ in cells]
                if min(xs) <= lo and max(xs) >= hi: continue
                for x, y in cells:
                    if lo <= x <= hi and (y < outer if above else y > outer):
                        outside += 1
                        break
        text = grid.render()
        lines = text.splitlines()
        # Feeds (or products) on rows right next to each other
        crowded = 0
        for drawn in labels_drawn.values():
            for i, (y1, a1, b1) in enumerate(drawn):
                for y2, a2, b2 in drawn[i + 1:]:
                    if abs(y1 - y2) == 1 and a1 <= b2 and a2 <= b1: crowded += 1 # Directly above or below
        # Junctions must have a line between them
        touching = 0
        chars = grid.chars
        for y in range(grid.height):
            row = chars[y]
            below = chars[y + 1] if y + 1 < grid.height else None
            for x in range(grid.width):
                if row[x] != '+': continue
                if x + 1 < width and row[x + 1] == '+': touching += 1
                if below is not None and below[x] == '+': touching += 1
        errors = overlaps + jogs + adjacent + outside + crowded + touching + grid.collisions + 10 * failed
        height = len(lines)
        width = max([len(i) for i in lines], default=0)
        return text, (errors, height, crossings, width, bends)

    def render(self, signs, extra):
        self.place_rows(signs)
        self.assign_ports()
        return self.draw(extra)

    def optimize(self):
        """
        Search for the best layout by iteratively moving to the best
        neighboring layout (flipping a loop or branch to the other side of
        the main line, or changing the spacing) until a layout repeats. The
        best layout found is the one with the fewest errors (overlapping
        lines, jogs, lines right next to each other, lines outside recycle
        loops they do not enclose, failed connections, and label
        collisions), then the smallest
        height, then the fewest crossings, the smallest width, and the
        fewest bends.
        """
        n_loops = self.n_loops
        default = tuple([-1] * n_loops + [1] * (len(self.groups) - n_loops))
        cache = {}
        def evaluate(state):
            if state not in cache:
                signs, spacing = state
                text, score = self.render(signs, spacing)
                deviations = sum([i != j for i, j in zip(signs, default)])
                cache[state] = (text, (*score, deviations))
            return cache[state][1]
        def neighbors(state):
            signs, spacing = state
            for i in range(len(signs)):
                yield (signs[:i] + (-signs[i],) + signs[i+1:], spacing)
            for other in SPACINGS:
                if other != spacing: yield (signs, other)
        state = (default, SPACINGS[0])
        visited = set()
        self.tight = set()
        self.row_extra = {}
        while state not in visited:
            visited.add(state)
            evaluate(state)
            state = min(neighbors(state), key=evaluate)
        best_state = min(cache, key=lambda s: cache[s][1])
        best_text, best_score = cache[best_state]
        # Pull units off the main line as close to it as possible, one group at a
        # time, as long as there are no more errors and the height does not grow
        signs, spacing = best_state
        # Only keep extra rows between rows of units where they are needed
        if spacing[0]:
            rows = sorted(set([n.row for n in self.nodes]))
            for r in rows[1:]:
                self.row_extra[r] = 0
                text, score = self.render(signs, spacing)
                if score[0] <= best_score[0]:
                    best_text, best_score = text, score
                else:
                    del self.row_extra[r]
        self.render(signs, spacing)
        for key in list(self.decoupled):
            self.tight.add(key)
            text, score = self.render(signs, spacing)
            if score[:2] <= best_score[:2]:
                best_text, best_score = text, score
            else:
                self.tight.discard(key)
        self.tight = set()
        self.row_extra = {}
        return best_text


def _stream_labels(IDs, length=5):
    """
    Return a dictionary of labels for feeds and products by stream ID: the
    first letters (5 by default) followed by any digits. Repeated labels are
    distinguished with numbers.
    """
    short = {}
    for ID in IDs:
        if not ID or ID in short: continue
        letters = ''.join([i for i in ID if i.isalpha()])[:length]
        digits = ''.join([i for i in ID if i.isdigit()])
        short[ID] = (letters, digits)
    groups = {}
    for ID, (letters, digits) in short.items():
        groups.setdefault(letters + digits, []).append(ID)
    labels = {}
    used = set([i for i, j in groups.items() if len(j) == 1])
    for label, group in groups.items():
        if len(group) == 1:
            labels[group[0]] = label
            continue
        number = 1
        for ID in group:
            letters = short[ID][0]
            while letters + str(number) in used: number += 1
            labels[ID] = letters + str(number)
            used.add(labels[ID])
            number += 1
    return labels


def ascii_flowsheet(nodes, edges, title=None, stream_label_length=5):
    """
    Return an ASCII diagram of a flowsheet.

    Parameters
    ----------
    nodes : Iterable[tuple[Hashable, str]]
        Key and label of each node. Node order should follow the
        processing order (edges going backwards are drawn as recycles).
    edges : Iterable[tuple[Hashable|None, int, str, int, Hashable|None]]
        Source key, source port, label, sink port, and sink key of each edge.
        A source of None represents a feed and a sink of None represents
        a product. Labels are only shown for feeds and products.
    title : str, optional
        Title placed above the diagram.
    stream_label_length : int, optional
        Maximum number of letters in feed and product labels (labels are the
        first letters of the stream name followed by any digits). Defaults to 5.

    """
    node_dct = {}
    for i, (key, label) in enumerate(nodes):
        node_dct[key] = _Node(key, label, i)
    node_list = list(node_dct.values())
    if not node_list: return title or ''
    edge_list = []
    pairs = set()
    edges = list(edges)
    labels = _stream_labels([e[2] for e in edges if (e[0] is None) != (e[4] is None)],
                            stream_label_length)
    edges = [(src, si, labels.get(label, label) if (src is None) != (dst is None) else label, di, dst)
             for src, si, label, di, dst in edges]
    for src, src_index, label, dst_index, dst in edges:
        src = node_dct.get(src)
        dst = node_dct.get(dst)
        if src is None and dst is None: continue
        label = label or ''
        if src is not None and dst is not None:
            if (src, dst) in pairs: continue # Parallel streams are drawn as one line
            pairs.add((src, dst))
        e = _Edge(src, src_index, label, dst_index, dst)
        if src is not None: src.outlets.append(e)
        if dst is not None: dst.inlets.append(e)
        edge_list.append(e)
    for n in node_list:
        n.inlets.sort(key=lambda e: e.dst_index)
        n.outlets.sort(key=lambda e: e.src_index)
        n.finalize()
    # Unconnected sections are drawn separately (units without connections
    # to other units are drawn together)
    sections = []
    singles = []
    seen = set()
    for n in node_list:
        if n in seen: continue
        section = []
        stack = [n]
        seen.add(n)
        while stack:
            m = stack.pop()
            section.append(m)
            for e in m.internal_inlets() + m.internal_outlets():
                other = e.src if e.dst is m else e.dst
                if other not in seen:
                    seen.add(other)
                    stack.append(other)
        if len(section) == 1:
            singles.extend(section)
        else:
            section.sort(key=lambda n: n.order)
            sections.append(section)
    if singles: sections.append(singles)
    diagrams = []
    for section in sections:
        members = set(section)
        section_edges = [e for e in edge_list if e.src in members or e.dst in members]
        diagrams.append(_Layout(section, section_edges).optimize())
    diagram = _pack(diagrams)
    if title: diagram = title + '\n\n' + diagram
    return diagram


def _pack(diagrams):
    """
    Place the largest diagram first and fit the others in its empty space
    (with a blank row above and below and two spaces on each side); diagrams
    that do not fit go below.
    """
    blocks = [i.split('\n') for i in diagrams if i]
    if not blocks: return ''
    blocks.sort(key=lambda b: -len(b) * max([len(i) for i in b]))
    canvas = [list(i) for i in blocks[0]]
    for block in blocks[1:]:
        H = len(canvas)
        W = max([len(i) for i in canvas])
        for row in canvas: row.extend([' '] * (W - len(row)))
        h = len(block)
        w = max([len(i) for i in block])
        # Prefix sums of occupied cells
        S = [[0] * (W + 1) for _ in range(H + 1)]
        for y in range(H):
            row = canvas[y]
            Sy = S[y + 1]
            Sp = S[y]
            total = 0
            for x in range(W):
                total += row[x] != ' '
                Sy[x + 1] = Sp[x + 1] + total
        def empty(y0, x0, y1, x1): # Inclusive bounds clipped to the canvas
            y0 = max(y0, 0); x0 = max(x0, 0); y1 = min(y1, H - 1); x1 = min(x1, W - 1)
            if y0 > y1 or x0 > x1: return True
            return S[y1 + 1][x1 + 1] - S[y0][x1 + 1] - S[y1 + 1][x0] + S[y0][x0] == 0
        position = None
        for y in range(H - h + 1):
            for x in range(W - w + 1):
                if empty(y - 1, x - 2, y + h, x + w + 1):
                    position = (y, x)
                    break
            if position: break
        if position:
            y, x = position
            for i, line in enumerate(block):
                for j, c in enumerate(line):
                    if c != ' ': canvas[y + i][x + j] = c
        else:
            canvas.append([])
            canvas.extend([list(i) for i in block])
    return '\n'.join([''.join(i).rstrip() for i in canvas])


# %% Interface with biosteam objects (biosteam is only imported when needed,
#    so that the layout engine above can be used without biosteam)

def _stream_label(stream):
    if not stream: return ''
    return stream.ID or ''

def _unit_label(unit, length=2):
    return unit.ID[:length].upper()

def _edges_from_connections(connections, units):
    import biosteam as bst
    show_all = bst.preferences.show_all_streams
    edges = []
    for source, source_index, stream, sink_index, sink in connections:
        has_source = source in units
        has_sink = sink in units
        if not (has_source or has_sink): continue
        if (has_source != has_sink and not show_all
            and stream and stream.isempty()): continue
        edges.append((
            source if has_source else None, source_index,
            _stream_label(stream),
            sink_index, sink if has_sink else None,
        ))
    return edges

def ascii_from_units(units, streams=None, title=None, unit_label_length=2, stream_label_length=5):
    """Return an ASCII diagram of the given units."""
    from biosteam.digraph.digraph import get_all_connections
    from biosteam.utils.piping import ignore_docking_warnings
    units = list(dict.fromkeys(units))
    if streams is None:
        streams = []
        stream_set = set()
        for u in units:
            for s in u._ins + u._outs:
                if s in stream_set: continue
                streams.append(s)
                stream_set.add(s)
    unit_set = set(units)
    connections = ignore_docking_warnings(get_all_connections)(streams)
    nodes = [(u, _unit_label(u, unit_label_length)) for u in units]
    edges = _edges_from_connections(connections, unit_set)
    return ascii_flowsheet(nodes, edges, title, stream_label_length)

def ascii_surface(path, title=None, unit_label_length=2, stream_label_length=5):
    """Return an ASCII diagram where subsystems are drawn as single units."""
    import biosteam as bst
    from biosteam.digraph.digraph import extend_surface_units, has_path
    from biosteam.utils.piping import ignore_docking_warnings
    surface_units = []
    old_unit_connections = set()
    isa = isinstance
    Unit = bst.Unit
    done = set()
    @ignore_docking_warnings
    def extend(i):
        extend_surface_units(i.ID, i.streams, i.units,
                             surface_units, old_unit_connections)
    for i in path:
        if i in done: continue
        done.add(i)
        if isa(i, Unit):
            surface_units.append(i)
        elif has_path(i):
            extend(i)
    try:
        return ascii_from_units(surface_units, title=title,
                                unit_label_length=unit_label_length,
                                stream_label_length=stream_label_length)
    finally:
        for u, ins, outs in old_unit_connections:
            u._ins[:] = ins
            u._outs[:] = outs

def ascii_minimal(ID, units, streams, title=None, stream_label_length=5):
    """Return an ASCII diagram of the system as a single unit."""
    import biosteam as bst
    from biosteam.digraph.digraph import get_section_inlets_and_outlets
    ins, outs = get_section_inlets_and_outlets(set(units), streams)
    show_all = bst.preferences.show_all_streams
    if not show_all:
        ins = [i for i in ins if not i.isempty()]
        outs = [i for i in outs if not i.isempty()]
    key = object()
    nodes = [(key, (ID or 'SYS').upper())]
    edges = [(None, 0, _stream_label(s), i, key) for i, s in enumerate(ins)]
    edges += [(key, i, _stream_label(s), 0, None) for i, s in enumerate(outs)]
    return ascii_flowsheet(nodes, edges, title, stream_label_length)
