from itertools import permutations
from multiprocessing import Pool, cpu_count

from dsenum.permutation_group import DerivativeStructurePermutation
from dsenum.coloring_generator import ColoringGenerator
from dsenum.core import hash_in_all_configuration, act_permutation


class DirectColoringEnumerator:
    """
    Parameters
    ----------
    permutation_group: list of permutation
    num_color: int
    cl_generator: ColoringGenerator
    color_exchange: bool
    """
    def __init__(self, permutation_group, num_color, cl_generator: ColoringGenerator,
                 color_exchange=True):
        self.permutation_group = permutation_group
        self.num_color = num_color
        self.cl_generator = cl_generator
        self.color_exchange = color_exchange

        self.list_colorings, self.flags = self.cl_generator.generate_all_colorings()

    def _hash(self, coloring):
        return hash_in_all_configuration(coloring, self.num_color)

    def _walk_orbit(self, coloring, include_identity=False):
        offset = 0 if include_identity else 1
        # assume self.permutation_group[0] is identity
        for prm in self.permutation_group[offset:]:
            acted_cl = act_permutation(prm, coloring)
            acted_cl_hash = self._hash(acted_cl)
            self.flags[acted_cl_hash] = False

    def coset_enumerate(self):
        colorings = []

        for cl in self.list_colorings:
            cl_hash = self._hash(cl)
            # avoid already-visited coloring
            if not self.flags[cl_hash]:
                continue
            colorings.append(cl)

            self._walk_orbit(cl, include_identity=False)

            if self.color_exchange:
                for cl_prm in permutations(range(self.num_color)):
                    exchanged_cl = [cl_prm[c] for c in cl]
                    exchanged_cl_hash = self._hash(exchanged_cl)
                    if exchanged_cl_hash == cl_hash:
                        continue
                    self._walk_orbit(exchanged_cl, include_identity=True)

        return colorings


class LexicographicColoringEnumerator:
    """
    this algorithm takes the most lexicographically small colorings

    Parameters
    ----------
    permutation_group: list of permutation
    num_color: int
    cl_generator: ColoringGenerator
    color_exchange: bool
    n_jobs: int
        when n_jobs > 1, use multiprocessing but it seems slower than signle core......
    """
    def __init__(self, permutation_group, num_color, cl_generator: ColoringGenerator,
                 color_exchange=True, n_jobs=1):
        self.permutation_group = permutation_group
        self.num_color = num_color
        self.cl_generator = cl_generator
        self.color_exchange = color_exchange
        if n_jobs == -1:
            self.n_jobs = cpu_count()
        else:
            self.n_jobs = n_jobs

    def _hash(self, coloring):
        return hash_in_all_configuration(coloring, self.num_color)

    def _is_champion_coloring(self, coloring):
        cl_hash = self._hash(coloring)

        if self.color_exchange:
            for cl_prm in permutations(range(self.num_color)):
                exchanged_cl = [cl_prm[c] for c in coloring]
                for prm in self.permutation_group:
                    acted_cl = act_permutation(prm, exchanged_cl)
                    acted_cl_hash = self._hash(acted_cl)
                    if acted_cl_hash < cl_hash:
                        return False
        else:
            for prm in self.permutation_group:
                acted_cl = act_permutation(prm, coloring)
                acted_cl_hash = self._hash(acted_cl)
                if acted_cl_hash < cl_hash:
                    return False
        return True

    def _is_champion_coloring_parallel(self, coloring):
        return coloring, self._is_champion_coloring(coloring)

    def coset_enumerate(self):
        if self.n_jobs != 1:
            with Pool(self.n_jobs) as pool:
                colorings = [cl for cl, keep
                             in pool.map(self._is_champion_coloring_parallel,
                                         self.cl_generator.yield_coloring())
                             if keep]
        else:
            colorings = []
            for cl in self.cl_generator.yield_coloring():
                if self._is_champion_coloring(cl):
                    colorings.append(cl)
        return colorings


class SiteColoringEnumerator(object):
    """
    Parameters
    ----------
    ds_permutation: DerivativeStructurePermutation object
    num_color: int
    color_exchange: bool
    leave_superperiodic: bool
    method: "direct" or "lexicographic", so far
    n_jobs: int, use only when method is "lexicographic"
    """

    def __init__(self, num_color,
                 ds_permutation: DerivativeStructurePermutation,
                 cl_generator: ColoringGenerator,
                 color_exchange=True, leave_superperiodic=False, use_all_colors=True,
                 method='direct', n_jobs=1):
        self.num_color = num_color
        self.ds_permutation = ds_permutation
        self.cl_generator = cl_generator
        self.color_exchange = color_exchange
        self.leave_superperiodic = leave_superperiodic
        self.use_all_colors = use_all_colors
        self.method = method
        self.n_jobs = n_jobs

        """
        if self.ds_permutation.num_site < self.num_color:
            raise ValueError('too many num_types')
        """

        self.permutation_group = self.ds_permutation.get_symmetry_operation_permutations()

        if self.method == 'direct':
            self.clenum = DirectColoringEnumerator(self.permutation_group, self.num_color,
                                                   self.cl_generator, color_exchange=color_exchange)
        elif self.method == 'lexicographic':
            self.clenum = LexicographicColoringEnumerator(self.permutation_group, self.num_color,
                                                          self.cl_generator,
                                                          color_exchange=color_exchange,
                                                          n_jobs=self.n_jobs)
        else:
            raise ValueError("Unknown method: ", self.method)

    @property
    def translation_permutations(self):
        return self.ds_permutation.prm_t

    def _hash(self, coloring):
        return hash_in_all_configuration(coloring, self.num_color)

    def unique_colorings(self):
        symmetry_uniqued_coloring = self.clenum.coset_enumerate()
        colorings = []

        for cl in symmetry_uniqued_coloring:
            if self.use_all_colors and (not self._has_all_colors(cl)):
                continue
            if (not self.leave_superperiodic) and self._is_superperiodic(cl):
                continue
            colorings.append(cl)

        return colorings

    def _has_all_colors(self, coloring):
        if len(set(coloring)) == self.num_color:
            return True
        else:
            return False

    def _is_superperiodic(self, coloring):
        # assume self.translation_permutations[0] is identity
        if any([self._hash(act_permutation(prm, coloring)) == self._hash(coloring)
               for prm in self.translation_permutations[1:]]):
            return True
        else:
            return False
