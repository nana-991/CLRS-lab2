from collections import deque
import sys

FORWARD = "FWD"
REVERSE = "REV"
CONTROL = "CTRL"

class SegmentNode:
    def __init__(self, kind, ref_idx, qry_idx):
        self.kind = kind
        self.ref_idx = ref_idx
        self.qry_idx = qry_idx

    def __hash__(self):
        return hash((self.kind, self.ref_idx, self.qry_idx))

    def __eq__(self, other):
        return (self.kind, self.ref_idx, self.qry_idx) == (other.kind, other.ref_idx, other.qry_idx)

class PathInfo:
    def __init__(self, cost, parent):
        self.cost = cost
        self.parent = parent

class GenomeMatcher:
    MAX_JUMP = 1000
    MIN_SEG_LEN = 15

    def __init__(self, query, reference):
        self.query = list(query)
        self.reference = list(reference)
        self.ref_rev = [self._comp(nuc) for nuc in reversed(reference)]

        self.memo = {}
        self.processed = set()
        self.queue = deque()

    def _comp(self, base):
        return {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}.get(base, base)

    def align(self):
        self._initialize()
        self._search_paths()
        self._print_results()

    def _initialize(self):
        start = SegmentNode(CONTROL, 0, 0)
        self.memo[start] = PathInfo(0, start)
        self.queue.append(start)

    def _search_paths(self):
        end_pos = len(self.query)
        while self.queue:
            current = self.queue.popleft()
            if current.qry_idx == end_pos:
                self.queue.appendleft(current)
                break
            if current in self.processed:
                continue
            self.processed.add(current)
            dist = self.memo[current].cost
            if current.kind == CONTROL:
                self._expand_control(current, dist)
            elif current.kind == FORWARD:
                self._expand_forward(current, dist)
            elif current.kind == REVERSE:
                self._expand_reverse(current, dist)

    def _expand_control(self, node, dist):
        q = node.qry_idx
        start = max(0, q - self.MAX_JUMP)
        end = min(len(self.reference), q + self.MAX_JUMP)
        for r in range(start, end + 1):
            self._try_update(SegmentNode(FORWARD, r, q), dist + 1, node)
            self._try_update(SegmentNode(REVERSE, r, q), dist + 1, node)

    def _expand_forward(self, node, dist):
        r, q = node.ref_idx, node.qry_idx
        if r < len(self.reference) and q < len(self.query):
            penalty = 0 if self.reference[r] == self.query[q] else 1
            self._try_update(SegmentNode(FORWARD, r + 1, q + 1), dist + penalty, node, penalty == 0)
        if r < len(self.reference):
            self._try_update(SegmentNode(FORWARD, r + 1, q), dist + 1, node)
        if q < len(self.query):
            self._try_update(SegmentNode(FORWARD, r, q + 1), dist + 1, node)
        self._try_update(SegmentNode(CONTROL, 0, q), dist + 1, node)

    def _expand_reverse(self, node, dist):
        r, q = node.ref_idx, node.qry_idx
        if r > 0 and q < len(self.query):
            penalty = 0 if self.ref_rev[r - 1] == self.query[q] else 1
            self._try_update(SegmentNode(REVERSE, r - 1, q + 1), dist + penalty, node, penalty == 0)
        if r > 0:
            self._try_update(SegmentNode(REVERSE, r - 1, q), dist + 1, node)
        if q < len(self.query):
            self._try_update(SegmentNode(REVERSE, r, q + 1), dist + 1, node)
        self._try_update(SegmentNode(CONTROL, 0, q), dist + 1, node)

    def _try_update(self, node, new_cost, parent, urgent=False):
        if node not in self.memo or self.memo[node].cost > new_cost:
            self.memo[node] = PathInfo(new_cost, parent)
            if urgent:
                self.queue.appendleft(node)
            else:
                self.queue.append(node)

    def _print_results(self):
        if not self.queue:
            print("[]")
            return
        node = self.queue[0]
        segs = []
        q_end = node.qry_idx
        r_end = node.ref_idx
        while not (node.kind == CONTROL and node.ref_idx == 0 and node.qry_idx == 0):
            prev = self.memo[node]
            if node.kind in {FORWARD, REVERSE} and node.kind != prev.parent.kind:
                q_start = node.qry_idx
                if node.kind == FORWARD:
                    r_start = node.ref_idx
                else:
                    r_start = r_end - 1
                if q_end - q_start + 1 > self.MIN_SEG_LEN:
                    segs.append((q_start, q_end, r_start, r_end))
            if node.kind == CONTROL and node.kind != prev.parent.kind:
                r_end = prev.parent.ref_idx
                q_end = prev.parent.qry_idx
            node = prev.parent
        print("[", end="")
        for i in range(len(segs) - 1, -1, -1):
            print(f"({segs[i][0]},{segs[i][1]},{segs[i][2]},{segs[i][3]})", end=", " if i > 0 else "")
        print("]")

if __name__ == "__main__":
    try:
        query_input = input().strip()
        reference_input = input().strip()
        matcher = GenomeMatcher(query_input, reference_input)
        matcher.align()
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
