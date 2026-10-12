"""FM-index: suffix array (SA-IS), Burrows-Wheeler transform, backward search.

Pure Python, zero dependencies. See README.md for the design notes.
"""
import argparse
import sys
from bisect import bisect_left


# ---------------------------------------------------------------- SA-IS
def sais(s, k):
    """Suffix array of int list `s` (alphabet 0..k-1) whose last symbol is a
    unique sentinel 0. Induced sorting, O(n)."""
    n = len(s)
    if n == 1:
        return [0]
    t = [False] * n  # True = S-type
    t[n - 1] = True
    for i in range(n - 2, -1, -1):
        t[i] = s[i] < s[i + 1] or (s[i] == s[i + 1] and t[i + 1])

    def is_lms(i):
        return i > 0 and t[i] and not t[i - 1]

    def buckets(ends):
        cnt = [0] * k
        for c in s:
            cnt[c] += 1
        out, tot = [], 0
        for c in cnt:
            tot += c
            out.append(tot if ends else tot - c)
        return out

    def induce(sa, lms_order):
        sa[:] = [-1] * n
        b = buckets(True)
        for i in reversed(lms_order):
            b[s[i]] -= 1
            sa[b[s[i]]] = i
        b = buckets(False)
        for j in range(n):
            i = sa[j] - 1
            if sa[j] > 0 and not t[i]:
                sa[b[s[i]]] = i
                b[s[i]] += 1
        b = buckets(True)
        for j in range(n - 1, -1, -1):
            i = sa[j] - 1
            if sa[j] > 0 and t[i]:
                b[s[i]] -= 1
                sa[b[s[i]]] = i

    lms = [i for i in range(1, n) if is_lms(i)]
    sa = [-1] * n
    induce(sa, lms)
    # name LMS substrings
    sorted_lms = [i for i in sa if is_lms(i)]
    names = [-1] * n
    name, prev = 0, -1
    for i in sorted_lms:
        diff = prev == -1
        if not diff:
            a, b = prev, i
            while True:
                if s[a] != s[b] or t[a] != t[b]:
                    diff = True
                    break
                a += 1
                b += 1
                if is_lms(a) or is_lms(b):
                    diff = not (is_lms(a) and is_lms(b) and s[a] == s[b])
                    break
        if diff:
            name += 1
        names[i] = name - 1
        prev = i
    reduced = [names[i] for i in lms]
    if name < len(lms):
        rsa = sais(reduced, name)
    else:
        rsa = [0] * len(lms)
        for i, c in enumerate(reduced):
            rsa[c] = i
    induce(sa, [lms[i] for i in rsa])
    return sa


def suffix_array(text):
    """Suffix array of a str/bytes (sentinel appended internally; result has
    len(text)+1 entries, first entry is len(text))."""
    if isinstance(text, str):
        text = text.encode("utf-8")
    syms = sorted(set(text))
    rank = {c: i + 1 for i, c in enumerate(syms)}
    s = [rank[c] for c in text] + [0]
    return sais(s, len(syms) + 1)


# ---------------------------------------------------------------- index
class FMIndex:
    """Full-text index over bytes. Query cost is independent of text length
    except for an O(log n) rank step."""

    def __init__(self, text, sa_step=16):
        if isinstance(text, str):
            text = text.encode("utf-8")
        if sa_step < 1:
            raise ValueError("sa_step must be >= 1")
        self.n = len(text)
        sa = suffix_array(text)
        self.sa_step = sa_step
        # BWT with sentinel represented as -1
        self.bwt = [text[i - 1] if i else -1 for i in sa]
        self.sampled = {r: p for r, p in enumerate(sa) if p % sa_step == 0}
        self._finish()

    def _finish(self):
        counts = {}
        self.pos = {}
        for i, c in enumerate(self.bwt):
            if c >= 0:
                counts[c] = counts.get(c, 0) + 1
                self.pos.setdefault(c, []).append(i)
        self.C, tot = {}, 1  # row 0 is the sentinel suffix
        for c in sorted(counts):
            self.C[c] = tot
            tot += counts[c]
        self.rows = self.n + 1

    def occ(self, c, i):
        """Occurrences of byte c in bwt[:i]."""
        p = self.pos.get(c)
        return bisect_left(p, i) if p else 0

    def _range(self, pattern):
        lo, hi = 0, self.rows
        for c in reversed(pattern):
            if c not in self.C:
                return 0, 0
            lo = self.C[c] + self.occ(c, lo)
            hi = self.C[c] + self.occ(c, hi)
            if lo >= hi:
                return 0, 0
        return lo, hi

    @staticmethod
    def _b(p):
        return p.encode("utf-8") if isinstance(p, str) else bytes(p)

    def count(self, pattern):
        lo, hi = self._range(self._b(pattern))
        return hi - lo

    def _resolve(self, row):
        steps = 0
        while row not in self.sampled:
            c = self.bwt[row]
            if c < 0:
                return steps
            row = self.C[c] + self.occ(c, row)
            steps += 1
        return self.sampled[row] + steps

    def locate(self, pattern):
        lo, hi = self._range(self._b(pattern))
        return sorted(self._resolve(r) for r in range(lo, hi))

    def extract(self, start, length):
        """Recover text[start:start+length] from the index alone."""
        end = min(self.n, start + length)
        if start < 0 or start > end:
            raise IndexError("bad range")
        # find row whose suffix begins at `end` by walking from a sample
        row = next((r for r, p in sorted(self.sampled.items(), key=lambda x: x[1])
                    if p >= end), None)
        out = []
        if row is None:
            row, p = 0, self.n  # sentinel suffix
        else:
            p = self.sampled[row]
        while p > start:
            c = self.bwt[row]
            out.append(c)
            row = self.C[c] + self.occ(c, row)
            p -= 1
        out.reverse()
        data = bytes(out)
        return data[: end - start]

    def invert(self):
        """Inverse BWT: rebuild the whole text (LF-walk from the sentinel)."""
        out, row = [], 0
        for _ in range(self.n):
            c = self.bwt[row]
            out.append(c)
            row = self.C[c] + self.occ(c, row)
        return bytes(reversed(out))

    def approx(self, pattern, k):
        """All (position, edit_distance_mismatches) with <= k substitutions,
        via backtracking backward search."""
        pat = self._b(pattern)
        res = {}
        alphabet = sorted(self.C)

        def rec(i, lo, hi, errs):
            if i < 0:
                for r in range(lo, hi):
                    p = self._resolve(r)
                    res[p] = min(errs, res.get(p, errs))
                return
            for c in alphabet:
                e = errs + (c != pat[i])
                if e > k:
                    continue
                nlo = self.C[c] + self.occ(c, lo)
                nhi = self.C[c] + self.occ(c, hi)
                if nlo < nhi:
                    rec(i - 1, nlo, nhi, e)

        rec(len(pat) - 1, 0, self.rows, 0)
        return sorted(res.items())

    # simple serialisation
    def dumps(self):
        import json
        return json.dumps({"n": self.n, "step": self.sa_step, "bwt": self.bwt,
                           "sa": sorted(self.sampled.items())})

    @classmethod
    def loads(cls, blob):
        import json
        d = json.loads(blob)
        self = cls.__new__(cls)
        self.n, self.sa_step, self.bwt = d["n"], d["step"], d["bwt"]
        self.sampled = {r: p for r, p in d["sa"]}
        self._finish()
        return self


def main(argv=None):
    ap = argparse.ArgumentParser(description="FM-index search")
    ap.add_argument("file")
    ap.add_argument("pattern")
    ap.add_argument("-k", type=int, default=0, help="allowed mismatches")
    ap.add_argument("--locate", action="store_true")
    a = ap.parse_args(argv)
    with open(a.file, "rb") as f:
        idx = FMIndex(f.read())
    if a.k:
        hits = idx.approx(a.pattern, a.k)
        print(len(hits), "hits")
        for p, e in hits[:50]:
            print(p, "mismatches=%d" % e)
    elif a.locate:
        for p in idx.locate(a.pattern):
            print(p)
    else:
        print(idx.count(a.pattern))
    return 0


if __name__ == "__main__":
    sys.exit(main())
