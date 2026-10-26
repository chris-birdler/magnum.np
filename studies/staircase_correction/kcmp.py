#!/usr/bin/env python
"""PLAN 16 E1: the compiled kernels of two inductor caches that started as copies of one base cache. For each file
(.py, .cpp, .ptx, .cubin, .so) by relative path: equal after the replacement of the cache path and the tree path by
placeholders (text files), or bytes equal (binary files; for a cubin also equal after the same replacement of the
path bytes). Files that exist in one cache only are listed with their normalized content hash, so that a kernel with
another name but the same content is matched.
    python kcmp.py CACHE_A TREE_A CACHE_B TREE_B"""
import hashlib, pathlib, sys
ca, ta, cb, tb = sys.argv[1:5]
EXT = {".py", ".cpp", ".ptx", ".cubin", ".so"}
def norm(p, cache, tree):
    b = pathlib.Path(p).read_bytes()
    return b.replace(cache.encode(), b"IC").replace((tree + "/").encode(), b"TREE/")
def files(c):
    return {str(f.relative_to(c)): f for f in pathlib.Path(c).rglob("*") if f.is_file() and f.suffix in EXT}
A, B = files(ca), files(cb)
same, diff, only = 0, [], []
for r in sorted(set(A) & set(B)):
    if norm(A[r], ca, ta) == norm(B[r], cb, tb):
        same += 1
    else:
        diff.append(r)
ha = {hashlib.md5(norm(A[r], ca, ta)).hexdigest(): r for r in set(A) - set(B)}
hb = {hashlib.md5(norm(B[r], cb, tb)).hexdigest(): r for r in set(B) - set(A)}
matched = set(ha) & set(hb)
print("%s / %s: %d / %d files; same relative path and content %d; same path, other content %d %s; "
      "only in one cache %d / %d, of them matched by content %d" % (ca, cb, len(A), len(B), same, len(diff), diff[:6],
      len(ha), len(hb), len(matched)))
for h in sorted(set(ha) - matched)[:4]: print("   only A:", ha[h])
for h in sorted(set(hb) - matched)[:4]: print("   only B:", hb[h])
