Supercell Construction
======================

.. role:: raw-latex(raw)
   :format: latex
..


Constructing supercell for a general transformation matrix
----------------------------------------------------------

We describe a rigorous way to associate a supercell with its primitive cell.
In what follows, we explain in a three-dimensional case.

We put a collection of point coordinates in a primitive cell as :math:`D`.
Any point coordinates :math:`\mathbf{x}` in the crystallographic pattern are uniquely decomposed as

.. math::
    :label: dsite

    \mathbf{x} = \mathbf{d} + \mathbf{m},\quad \mathbf{d} \in D,\, \mathbf{m} \in \mathbb{Z}^{3}

where :math:`\mathbf{m}` indicates indices of a lattice point of a unit cell that the point coordinates :math:`\mathbf{x}` belongs.
When we consider translational symmetry of a sublattice :math:`L_{\mathbf{M}}`, two lattice points :math:`\mathbf{m}, \mathbf{m}'` are equivalent if the distance between the two lattice points is a translation of :math:`L_{\mathbf{M}}`, that is,

.. math::
    :label: eq:equiv-lattice-points

    \mathbf{m} - \mathbf{m}' \in \mathbf{M}\mathbb{Z}^{3}.

The Smith normal form (SNF) of the transformation matrix :math:`\mathbf{M}` is useful to concretely write down Eq. :eq:`eq:equiv-lattice-points`.
The SNF is one of the decomposition of an integer matrix :math:`\mathbf{M}` as

.. math::
    :label: eq:SNF

    \mathbf{S} = \mathbf{PMQ}

where :math:`\mathbf{P}` and :math:`\mathbf{Q}` are unimodular matrices, and :math:`\mathbf{S}` is a diagonal integer matrix,

.. math::

    \mathbf{S} =
        \begin{pmatrix}
            S_{11} & 0 & 0 \\
            0 & S_{22} & 0 \\
            0 & 0 & S_{33}
        \end{pmatrix}.

Here :math:`S_{11}` is a divisor of :math:`S_{22}`, and :math:`S_{22}` is a divisor of :math:`S_{33}`.
We can rewrite Eq. :eq:`eq:equiv-lattice-points` with Eq. :eq:`eq:SNF` as

.. math::
    :label: eq:equiv-represents

    [\mathbf{Pm}]_{\mathbf{S}} = [\mathbf{Pm}']_{\mathbf{S}}

where :math:`[\cdot]_{\mathbf{S}}` indicates to take modulus for the :math:`i`\ th row by :math:`S_{ii}`.
We mention that the range of :math:`[\cdot]_{\mathbf{S}}` is :math:`\mathbb{Z}_{S_{11}} \oplus \mathbb{Z}_{S_{22}} \oplus \mathbb{Z}_{S_{33}}` because a value of the :math:`i`\ th row is a remainder by :math:`S_{ii}`.
Using Eq. :eq:`eq:equiv-represents`, the point coordinates :math:`\mathbf{x}` are decomposed with the translation of the sublattice :math:`L_{\mathbf{M}}` as

.. math::
    :label: eq:site-decomposition

    \mathbf{x} = \mathbf{d} + \mathbf{P}^{-1} \mathbf{f} + \mathbf{Ml} ,\, \quad \mathbf{f} \in \mathbb{Z}_{S_{11}} \oplus \mathbb{Z}_{S_{22}} \oplus \mathbb{Z}_{S_{33}},\, \mathbf{l} \in \mathbb{Z}^{3}.

In Eq. :eq:`eq:site-decomposition`, the terms :math:`\mathbf{d} + \mathbf{P}^{-1}\mathbf{f}` denote point coordinates in the supercell and :math:`\mathbf{l}` stands for the indices of a supercell.

The collection of point coordinates in :math:`D_{\mathbf{M}}` is the set of sites inside the parallelepiped spanned by the set of basis vectors :math:`\mathbf{AM}`,

.. math::

    D_{\mathbf{M}} = \left\{ \mathbf{d} + \mathbf{m} \mid \mathbf{d} \in D, \mathbf{m} \in \mathbb{Z}^{3} \right\} \cap \mathbf{M} [0, 1)^{3} .

Using Eq. :eq:`eq:site-decomposition`, we have explicit point coordinates in the supercell:

.. math::
    :label: eq:concrete-supercell

    D_{\mathbf{M}} =
    \left\{
        \mathbf{d} + \mathbf{P}^{-1} \mathbf{f} + \mathbf{M} \mathbf{l}_{\mathbf{f}} \mid
        \mathbf{d} \in D,\, \mathbf{f} \in \mathbb{Z}_{S_{11}} \oplus \mathbb{Z}_{S_{22}} \oplus \mathbb{Z}_{S_{33}}
    \right\}

where :math:`\mathbf{l}_{\mathbf{f}}` is an offset to locate each site in :math:`D_{\mathbf{M}}` inside the parallelepiped.

Then we define a mapping that embeds any point coordinates :math:`\mathbf{x} = \mathbf{d} + \mathbf{m} \, (\mathbf{d} \in D, \mathbf{m} \in \mathbb{Z}^{3})` in the supercell as

.. math::

    \pi_{\mathbf{M}}(\mathbf{d} + \mathbf{m}) = \mathbf{d} + \mathbf{P}^{-1} [\mathbf{Pm}]_{\mathbf{S}} + \mathbf{M} \mathbf{l}_{[\mathbf{Pm}]_{\mathbf{S}}} \in D_{\mathbf{M}}.

In what follows, we label the point coordinates in the supercell as :math:`D_{\mathbf{M}} = \left\{ \mathbf{d}_{1}, \dots, \mathbf{d}_{|D_{\mathbf{M}}|} \right\}`.

Putting a space group of the supercell :math:`(\mathbf{AM}, D_{\mathbf{M}})` as :math:`\mathcal{H}_{\mathbf{M}}` and a space group of its primitive cell as :math:`\mathcal{G}`, the space group :math:`\mathcal{H}_{\mathbf{M}}` is a subgroup of the space group :math:`\mathcal{G}`.
When :math:`g = \{ \mathbf{R} \mid \mathbf{\tau} \} \in \mathcal{G}` is contained in :math:`\mathcal{H}_{\mathbf{M}}`, it is required that the sublattice :math:`L_{\mathbf{M}}` and the rotated one :math:`L_{\mathbf{RM}}` are coincided, that is, :math:`\mathbf{M}^{-1} \mathbf{RM}` is unimodular.
By applying each symmetry operation :math:`{ g = \{ \mathbf{R} \mid \mathbf{\tau} \} \in \mathcal{H}_{\mathbf{M}} }`, point coordinates :math:`\mathbf{d}_{i}` are moved to :math:`g \mathbf{d}_{i} (= \mathbf{R}\mathbf{d}_{i} + \mathbf{\tau})` and the latter point coordinates are equivalent to :math:`\pi_{\mathbf{M}} (g \mathbf{d}_{i})` up to translations in the sublattice :math:`L_{\mathbf{M}}`.
Thus a permutation representation :math:`\sigma_{g}` of operation :math:`g` is obtained by

.. math::

    \mathbf{d}_{\sigma_{g}(i)} = \pi_{\mathbf{M}} (g \mathbf{d}_{i}).

Although :math:`\mathcal{H}_{\mathbf{M}}` is infinite group, its permutation group :math:`\Sigma_{\mathbf{M}} = \{ \sigma_{g} \mid g \in \mathcal{H}_{\mathbf{M}} \}` is finite.
Putting a translation subgroup obtained from a lattice :math:`L` as

.. math::

    T_{L} = \left\{ \{ \mathbf{I} \mid \mathbf{t} \} \mid \mathbf{t} \in L \right\} \subseteq \mathcal{H}_{\mathbf{M}},

we only have to calculate the permutation representation for its factor group :math:`T_{L} / T_{L_{\mathbf{M}}}` because translations from the sublattice :math:`L_{\mathbf{M}}` trivially return the identity permutation for a collection of point coordinates :math:`D_{\mathbf{M}}`.
We can construct one representative for :math:`T_{L} / T_{L_{\mathbf{M}}}` in the same way as Eq. :eq:`eq:concrete-supercell`,

.. math::

    T_{L} / T_{L_{\mathbf{M}}} =
    \left\{
        \{\mathbf{I} \mid \mathbf{P}^{-1} \mathbf{f} \} \mid
        \mathbf{f} \in \mathbb{Z}_{S_{11}} \oplus \mathbb{Z}_{S_{22}} \oplus \mathbb{Z}_{S_{33}}
    \right\}.

The permutation group :math:`\Sigma_{\mathbf{M}}` is obtained from finite operations in :math:`\mathcal{H}_{\mathbf{M}} / T_{L}` and :math:`T_{L} / T_{L_{\mathbf{M}}}` as its generators.

Convention in dsenum
--------------------
`dsenum.site.DerivativeSite` represents a pair of :math:`\mathbf{d}` and :math:`\mathbf{m}` in Eq. :eq:`dsite`.

The following instance of `DerivativeSite`

.. code-block:: python

    dsite = dsenum.site.DerivativeSite(site_index, jimage)

corresponds to

.. math::
    \mathbf{x} = \mathbf{d}_{\mathrm{site\_index}} + \mathrm{jimage}.

`dsenum.site.CanonicalSite` represents a pair of :math:`\mathbf{d}` and :math:`\mathbf{f}` in Eq. :eq:`eq:site-decomposition`.

The following instance of `CanonicalSite`

.. code-block:: python

    csite = dsenum.site.CannonicalSite(site_index, f)

corresponds to

.. math::
    \mathbf{x} = \mathbf{d}_{\mathrm{site\_index}} + \mathbf{P}^{-1} \mathbf{f} + \mathbf{Ml}

References
----------
* Gus L. W. Hart and Rodney W. Forcade, "Algorithm for generating derivative structures," Phys. Rev. B 77 224115, (2008)
* Gus L. W. Hart and Rodney W. Forcade, "Generating derivative structures from multilattices: Application to hcp alloys," Phys. Rev. B 80 014120 (2009)
* Lyuwen Fu, Mordechai Kornbluth, Zhengqian Cheng, and Chris A. Marianetti, "Group theoretical approach to computing phonons and their interactions", Phys. Rev. B 100, 014303 (2019)
