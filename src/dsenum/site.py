from dataclasses import dataclass


@dataclass(unsafe_hash=True)
class DerivativeSite:
    # jimage indicates index of unitcell
    site_index: int
    jimage: tuple[int, ...]


@dataclass(unsafe_hash=True)
class CanonicalSite:
    # factor indicates index of lattice points with some periodic condition
    site_index: int
    factor: tuple[int, ...]
