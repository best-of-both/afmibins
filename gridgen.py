#!python
from math import sqrt

def generate_grid(maxsize=0.1):
    phi = (1 + sqrt(5)) / 2;
    den = sqrt(sqrt(5) * phi);

    pts = [
        (0, -den, 0),
        (0,  phi / den, 2 * phi / den),
        (0, -phi / den, -2 * phi / den),
        (0, den, 0),
        (phi, -phi / den, -1 / den),
        (-phi, -phi / den, -1 / den),
        (phi, phi / den, 1 / den),
        (-phi, phi / den, 1 / den),
        (1, -phi / den, phi * phi / den),
        (1, phi / den, -phi * phi / den),
        (-1, -phi / den, phi * phi / den),
        (-1, phi / den, -phi * phi / den)
    ]

    tr = [
        (0, 2, 4),
        (0, 2, 5),
        (0, 4, 8),
        (0, 5, 10),
        (0, 8, 10),
        (1, 3, 6),
        (1, 3, 7),
        (1, 6, 8),
        (1, 7, 10),
        (1, 8, 10),
        (2, 4, 9),
        (2, 5, 11),
        (2, 9, 11),
        (3, 6, 9),
        (3, 7, 11),
        (3, 9, 11),
        (4, 6, 9),
        (4, 6, 8),
        (5, 7, 10),
        (5, 7, 11)
    ]

    size = den
    while size > den * maxsize:
        new_pts = {}
        new_size = size
        new_tr = []
        for a, b, c in tr:
            pa = pts[a]
            pb = pts[b]
            pc = pts[c]
            ab_i = get_or_create_pt(a, b, new_pts, pts, den)
            ac_i = get_or_create_pt(a, c, new_pts, pts, den)
            bc_i = get_or_create_pt(b, c, new_pts, pts, den)
            ab = pts[ab_i]
            ac = pts[ac_i]
            bc = pts[bc_i]
            new_tr.append((a, ab_i, ac_i))
            new_tr.append((b, ab_i, bc_i))
            new_tr.append((c, ac_i, bc_i))
            new_tr.append((ab_i, ac_i, bc_i))
            ps = max(dist(pa, ab), dist(pb, bc), dist(pc, ac),
                     dist(ab, bc), dist(ab, ac), dist(ac, bc))
            if ps < new_size:
                new_size = ps
        tr = new_tr
        size = new_size
    return pts, size / den

def dist(a, b):
    return sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2 + (a[2] - b[2])**2)

def get_or_create_pt(a, b, new_pts, pts, den):
    if (a, b) in new_pts:
        return new_pts[(a, b)]
    pa = pts[a]
    pb = pts[b]
    ab = ((pa[0] + pb[0]) / 2,
          (pa[1] + pb[1]) / 2,
          (pa[2] + pb[2]) / 2)
    ab_r = den / dist(ab, (0, 0, 0))
    ab = (ab_r * ab[0], ab_r * ab[1], ab_r * ab[2])
    ab_i = len(pts)
    new_pts[(a, b)] = ab_i
    pts.append(ab)
    return ab_i

if __name__ == '__main__':
    import sys
    maxsize = None
    if len(sys.argv) > 1:
        maxsize = float(sys.argv[1])
    grid, size = generate_grid(maxsize)
    print('[')
    for p in grid:
        print(p[0], p[1], p[2])
    print(']')
    print('%', len(grid), size)
