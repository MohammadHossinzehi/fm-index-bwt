# fm-index-bwt

A compact full-text search index in pure Python: **linear-time suffix array construction (SA-IS), the Burrows-Wheeler transform, and FM-index backward search** with sampled locate, text extraction, inverse BWT, and k-mismatch approximate matching. Zero dependencies.

## Why it is useful

A naive substring search scans the whole text per query. An FM-index (the structure behind read aligners like Bowtie and BWA) answers `count(pattern)` in O(m log n) time for a pattern of length m, regardless of text size, and can even reproduce the original text, so the index *replaces* the text.

## What is in the box

- `suffix_array(text)`: SA-IS induced sorting, O(n).
- `FMIndex(text, sa_step=16)`: builds the BWT, `C` table and per-symbol rank structure.
  - `count(p)`, `locate(p)`: exact matching by backward search; locate walks LF-mappings to a sampled suffix array entry (time/space trade-off via `sa_step`).
  - `extract(start, length)`: recover a substring from the index alone.
  - `invert()`: inverse BWT, rebuilds the entire text.
  - `approx(p, k)`: all positions matching with at most k substitutions, using backtracking over the backward search with pruning of empty ranges.
  - `dumps()` / `loads()`: JSON serialisation.
- A small CLI.

## Run it

```
python3 fmindex.py corpus.txt "needle"            # count
python3 fmindex.py corpus.txt "needle" --locate   # positions
python3 fmindex.py corpus.txt "needle" -k 1       # up to 1 mismatch
python3 -m unittest -v                            # tests
```

```python
from fmindex import FMIndex
idx = FMIndex("banana", sa_step=2)
idx.count("ana")      # 2
idx.locate("ana")     # [1, 3]
idx.invert()          # b'banana'
```

## Design decisions

- Texts are treated as bytes, so UTF-8 works naturally. A unique sentinel sorts before every symbol; internally it is symbol 0 for SA-IS and `-1` in the BWT.
- Rank (`occ`) uses a sorted position list per symbol plus `bisect`: O(log n) and very little code. A wavelet tree or checkpointed bit vectors would be the next step for memory.
- Suffix array sampling keeps every `sa_step`-th text position. Unsampled rows are resolved by LF-stepping, so locate costs O(sa_step) rank calls per hit.

## Testing

SA-IS is checked against a naive sort on hundreds of random strings and on highly repetitive inputs (`aaaa...`, periodic strings). Count, locate, extract, inverse BWT and approximate search are compared with brute force for several sampling rates, plus serialisation and Unicode cases.
