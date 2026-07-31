Definitions and Terms
=====================

This page is a work in progress.

The calculations implemented in `stellar_geology` are from Putirka & Rarick (2019) and benchmarked to the manuscript's supplementary data files.

Variable Definitions
--------------------
Acronyms used throughout this documentation:
    - BP = Bulk Planet
    - BSP = Bulk Silicate Planet (BP minus the core)
    - CMF = Core mass fraction

alphas
^^^^^^
    Alpha values as described in Putirka and Rarick (2019) represent the the proportion of an element retained within the mantle relative to the bulk planet. This value can be set for any element known to `stellar_geology` For example:

    .. math::

        \alpha{Fe} = Fe^{BSP}/Fe^{BP}

    where :math:`Fe^{BSP}` is Fe in a bulk silicate planet (bulk planet, minus core) on a cation weight % basis (elemental weight proportions, absent anions) and :math:`Fe^{BP}` is the cation weight % of Fe for the bulk planet.

    Putirka & Rarick (2019) assume that :math:`\alpha_{Fe}^{Mercury}` = 0.0–0.12; :math:`\alpha_{Fe}^{Earth}` = 0.263–0.494; and :math:`\alpha_{Fe}^{Mars}` = 0.54–0.58.

    An alpha is the ratio of an element's concentration in the BSP to its concentration in the bulk planet, with both compositions normalized to 100%. It is *not* the fraction of that element's mass residing in the mantle. Removing core mass concentrates everything left behind, so the two differ by the silicate mass fraction: mass retained in mantle = alpha × `silicate_mass_fraction`. For example, with alpha_Fe = 0.49 and a 10% core, only ~44% of the planet's Fe mass is in the mantle — not 49%. Quoting alphas as retention fractions overstates mantle retention. The lithophile enrichment itself is informative: the common BSP/BP ratio of the fully lithophile elements fixes the core size, which is how `Planet.core_mass_fraction` and `Planet.core_composition` are derived.

    For a model-consistent pair, all fully lithophile elements share one enrichment factor exactly. A measured or independently modeled BSP generally violates this, and then no alpha set can exactly reproduce the pair. In this case, the derived alphas and core properties are best-fit approximations, and `stellar_geology`` raises a `UserWarning` when the lithophile ratios disagree beyond a set tolerance.

    **On the reproducability of alphas**

    In `stellar_geology`, a BP + BSP pair can be used to calculate a set of alpha values. In theory, a perfect solution would mean that deriving a new BSP from the original BP and calculated alphas should exactly match the original BSP. In practice, there are several reasons why this is not true. For a BSP to be perfectly reproducable, the following must also be true:
    
    1. The BP and BSP must be perfectly matched. That is, there must be exactly one set of alpha values that mathematically connects the two. In forward modeling with `stellar_geology`, a BSP calculated from a BP and alphas is perfectly matched to that BP. However, if alpha values are calculated from a BP + BSP pair from two difference sources (e.g., measured stellar composition and published BSP composition), there may be no alpha set that satisfies the relationship exactly. If any mismatch is found during computation that is much larger than floating point error, `stellar_geology` will warn you that the core frationation parameters calculated, like alphas, are best fit approximations.
    2. The model must consider volatile elements. Currently, `stellar_geology` does not consider any volatile elements, the presence of which would variably dilute the refractory composition, since all compositions are normalized to 100%.
    3. The model must consider mineralogy. Currently, `stellar_geology` partitions elements between the mantle and core on a per-element basis. It does not consider oxygen or any other anion that would in reality be bound to those elements in the silicate. For example, Fe is metallic Fe in the core but bound in minerals like olivine and pyroxene in the mantle. The alpha value is a useful workaround for this for simple calculations like those performed in `stellar_geology` in that it is mineralogy agnostic, allowing us to make larger inferences about bulk compositions without the need for complicated phase modeling. Although this is a significant assumption we make, it strips away several complicating factors that may be unneccessary when examining planetary-scale processes, particularly on exoplanets.


core_mass_fraction
^^^^^^^^^^^^^^^^^^
    At present, the `core_mass_fraction` in `stellar_geology` is only a derived value and is computed from a BP + BSP pair, using the same math and assumptions for calculating alphas from a BP + BSP pair. `core_mass_fraction` is a mass fraction on the volatile-free, oxygen-free element basis; a true planetary CMF additionally requires the oxygen bound in mantle oxides and any light elements in the core.
    
    An alpha is the ratio of an element's concentration in the BSP to its concentration in the bulk planet, with both compositions normalized to 100%. It is *not* the fraction of that element's mass residing in the mantle. Removing core mass concentrates everything left behind, so the two differ by the silicate mass fraction: mass retained in mantle = alpha × `silicate_mass_fraction`. For example, with alpha_Fe = 0.49 and a 10% core, only ~44% of the planet's Fe mass is in the mantle — not 49%. Quoting alphas as retention fractions overstates mantle retention. The lithophile enrichment itself is informative: the common BSP/BP ratio of the fully lithophile elements fixes the core size, which is how `Planet.core_mass_fraction` and `Planet.core_composition` are derived.

    BSP and BP are normalized compositions, which only constrain alphas relative to the most silicate-enriched elements. Deriving alphas (and core properties) from a BP + BSP pair assumes those most-enriched elements are perfectly lithophile — a uniform depletion of every element into the core is invisible to normalized compositions. This is the same convention used by Putirka & Rarick (2019), but it is an assumption, not something the compositions can test.

References
----------
Putirka, K.D. and Rarick, J.C. (2019) The Composition and Mineralogy of Rocky Exoplanets: A Survey of >4,000 Stars from the Hypatia Catalog