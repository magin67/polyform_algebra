"""Exact continuants for open chains of cycle cores."""


def cycle_chain_tree_count(core_sizes, seam_lengths):
    """Return the spanning-tree count of an admissible open cycle chain.

    Sizes count unit edges on the cycles. Seam lengths count shared edges,
    so a one-vertex seam has length zero. The input describes sizes only.
    The actual assembly must have disjoint seams on each intermediate core.
    """
    sizes = tuple(core_sizes)
    seams = tuple(seam_lengths)
    if not sizes or any(type(size) is not int or size < 3 for size in sizes):
        raise ValueError("expected cycle sizes of at least three")
    if len(seams) != len(sizes) - 1 or any(type(length) is not int or length < 0 for length in seams):
        raise ValueError("expected one nonnegative integer length per seam")
    if any(left + right + int(i > 0) + int(i < len(sizes) - 1) > size
           for i, (size, left, right) in enumerate(zip(sizes, (0, *seams), (*seams, 0), strict=True))):
        raise ValueError("seam vertex counts exceed a core size")
    if any(length >= min(left, right) - 1
           for left, right, length in zip(sizes[:-1], sizes[1:], seams, strict=True)):
        raise ValueError("a seam must be a proper induced path in each core")
    previous, current = 1, sizes[0]
    for size, length in zip(sizes[1:], seams, strict=True):
        previous, current = current, size * current - length**2 * previous
    return current
