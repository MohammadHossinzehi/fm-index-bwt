import random
import unittest
from fmindex import FMIndex, suffix_array


def naive_sa(b):
    if isinstance(b, str):
        b = b.encode()
    return [len(b)] + sorted(range(len(b)), key=lambda i: b[i:])


def naive_find(t, p):
    t = t.encode() if isinstance(t, str) else t
    p = p.encode() if isinstance(p, str) else p
    return [i for i in range(len(t) - len(p) + 1) if t[i:i + len(p)] == p]


class T(unittest.TestCase):
    def test_banana(self):
        self.assertEqual(suffix_array("banana"), [6, 5, 3, 1, 0, 4, 2])
        idx = FMIndex("banana", sa_step=2)
        self.assertEqual(idx.count("ana"), 2)
        self.assertEqual(idx.locate("ana"), [1, 3])
        self.assertEqual(idx.count("x"), 0)

    def test_sa_random(self):
        rnd = random.Random(1)
        for _ in range(300):
            n = rnd.randint(1, 60)
            t = bytes(rnd.choice(b"ab" if _ % 2 else b"abcd") for _ in range(n))
            self.assertEqual(suffix_array(t), naive_sa(t))

    def test_repetitive(self):
        for t in [b"a" * 50, b"ab" * 40, b"abaababaabaab" * 5]:
            self.assertEqual(suffix_array(t), naive_sa(t))

    def test_queries_random(self):
        rnd = random.Random(2)
        for step in (1, 3, 16):
            t = bytes(rnd.choice(b"acgt") for _ in range(400))
            idx = FMIndex(t, sa_step=step)
            for _ in range(60):
                i = rnd.randrange(len(t))
                p = t[i:i + rnd.randint(1, 8)]
                self.assertEqual(idx.locate(p), naive_find(t, p))
                self.assertEqual(idx.count(p), len(naive_find(t, p)))

    def test_invert_and_extract(self):
        t = b"the quick brown fox jumps over the lazy dog" * 3
        idx = FMIndex(t, sa_step=5)
        self.assertEqual(idx.invert(), t)
        rnd = random.Random(3)
        for _ in range(50):
            s = rnd.randrange(len(t))
            l = rnd.randint(0, 30)
            self.assertEqual(idx.extract(s, l), t[s:s + l])

    def test_approx(self):
        rnd = random.Random(4)
        t = bytes(rnd.choice(b"acgt") for _ in range(300))
        idx = FMIndex(t, sa_step=4)
        for _ in range(20):
            i = rnd.randrange(len(t) - 10)
            p = bytearray(t[i:i + 10])
            p[rnd.randrange(10)] = ord("a")
            p = bytes(p)
            exp = {}
            for j in range(len(t) - 9):
                e = sum(x != y for x, y in zip(t[j:j + 10], p))
                if e <= 2:
                    exp[j] = e
            self.assertEqual(dict(idx.approx(p, 2)), exp)

    def test_roundtrip_serialisation(self):
        idx = FMIndex("mississippi", sa_step=2)
        j = FMIndex.loads(idx.dumps())
        self.assertEqual(j.locate("ssi"), [2, 5])

    def test_unicode_and_empty_hit(self):
        idx = FMIndex("héllo wörld héllo")
        self.assertEqual(idx.count("héllo"), 2)


if __name__ == "__main__":
    unittest.main()
