Definitions and Terms
=====================

This page is a work in progress.

The calculations implemented in the direct Star -> Planet pipeline in `stellar_geology` are from Putirka & Rarick (2019) and benchmarked to the manuscript's supplementary data files. We define additional equations to derive mantle enrichment of lithophile elements and core mass fraction.

Variable Definitions
--------------------
Acronyms used throughout this documentation:
    - BP = Bulk Planet
    - BSP = Bulk Silicate Planet (BP minus the core)
    - CMF = Core mass fraction
    - SMF = Silicate mass fraction (1-CMF)
  
Variable names:
    - :math:`m`: mass retention factor
    - :math:`\alpha`: "alpha", the concentration ratio of an element between BSP and BP
    - :math:`\phi`: the lithophile enrichment factor

Other terms:
    - lithophile: an element that partitions strongly into a silicate phase relative to a metal or vapor phase. Here, the term "lithophile" is used to describe elements where :math:`\alpha` > 1, meaning that they become enriched in the mantle after core formation.

alphas
^^^^^^
    Alpha values as described in Putirka and Rarick (2019) represent the the proportion of an element retained within the mantle relative to the bulk planet. This value can be set for any element known to `stellar_geology` For example:

    .. math::

        \alpha{Fe} = [Fe]^{BSP}/[Fe]^{BP}

    where :math:`[Fe]^{BSP}` is the concentration of Fe in a bulk silicate planet (bulk planet, minus core) on a cation weight % basis (elemental weight proportions, absent anions) and :math:`[Fe]^{BP}` is the cation weight % of Fe for the bulk planet. 

    Putirka & Rarick (2019) assume that :math:`\alpha_{Fe}^{Mercury}` = 0.0–0.12; :math:`\alpha_{Fe}^{Earth}` = 0.263–0.494; and :math:`\alpha_{Fe}^{Mars}` = 0.54–0.58.

    **Derived alphas**

    For any element that partitions into the core (is not fully lithophile), alpha are contained in the set (0,1]. Where :math:`\alpha` is derived from known BP and BSP compositions, :math:`\alpha` = 0 means that the element is missing from either the BSP or BP; it cannot be represented using the alpha variable. Derived :math:`\alpha` >1 means an element is strongly lithophile relative to elements with :math:`\alpha` ≤ 1. See the explanation of the lithophile enrichment factor below.

    **On the reproducability of alphas**

    In `stellar_geology`, a BP + BSP pair can be used to calculate a set of alpha values. In theory, a perfect solution would mean that a derived alpha and the original known BP would recover the original known BSP exactly. In practice, there are several reasons why this is not true. For a BSP to be perfectly reproducable, the following must also be true:
    
    1. The BP and BSP must be perfectly matched. That is, there must be exactly one set of alpha values that mathematically connects the two. In forward modeling with `stellar_geology`, a BSP calculated from a BP and alphas is perfectly matched to that BP. However, if alpha values are calculated from a BP + BSP pair from two difference sources (e.g., measured stellar composition and published BSP composition), there may be no alpha set that satisfies the relationship exactly. If any mismatch is found during computation that is much larger than floating point error, `stellar_geology` will warn you that the core frationation parameters calculated, like alphas, are best fit approximations.
    2. The model must consider volatile elements. Currently, `stellar_geology` does not consider any volatile elements, the presence of which would variably dilute the refractory composition, since all compositions are normalized to 100%.
    3. The model must consider mineralogy. Currently, `stellar_geology` partitions elements between the mantle and core on a per-element basis. It does not consider oxygen or any other anion that would in reality be bound to those elements in the silicate. For example, Fe is metallic Fe in the core but bound in minerals like olivine and pyroxene in the mantle. The alpha value is a useful workaround for this for simple calculations like those performed in `stellar_geology` in that it is mineralogy agnostic, allowing us to make larger inferences about bulk compositions without the need for complicated phase modeling. Although this is a significant assumption we make, it strips away several complicating factors that may be unneccessary when examining planetary-scale processes, particularly on exoplanets.

Lithophile enrichment factor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    An alpha is the ratio of an element's concentration in the BSP to its concentration in the bulk planet, with both compositions normalized to 100%. It is *not* equal to the mass ratio of that element in the BSP. Removing core mass concentrates everything left behind in the mantle, and so strongly lithophile elements can have :math:`\alpha` > 1. For such elements, we define a lithophile enrichment factor described by:

    .. math::

        \phi = \frac{\Sigma[i_{lithophile}]^{BSP}}{\Sigma[i_{lithophile}]^{BP}}
    
    where :math:`\phi` is the lithophile enrichment factor shared by all perfectly lithophile elements, and :math:`[i_{lithophile}]` is the concentration of component :math:`i` in the BSP or BP.

    For a model-consistent BP-BSP pair, all fully lithophile elements are enriched to the same degree as each other and so share a common :math:`\phi` value. The lithophile enrichment itself is informative: the common :math:`\phi` of the fully lithophile elements fixes the core size, which is how CMF (`Planet.core_mass_fraction`) and `Planet.core_composition` are derived (see below).
    
    A measured or independently modeled BSP generally violates this, and then no alpha set can exactly reproduce the pair. In this case, the derived alphas and core properties are best-fit approximations, and `stellar_geology`` raises a `UserWarning` when the lithophile ratios disagree beyond a set tolerance.

core_mass_fraction
^^^^^^^^^^^^^^^^^^
    At present, the `core_mass_fraction` in `stellar_geology` is only a derived value and is computed from a BP + BSP pair, using the same assumptions for calculating alphas from a BP + BSP pair. We can calculate the CMF and SMF directly from the lithophile enrichment factor :math:`\phi` as:

    .. math::

        CMF = \frac{\phi - 1}{\phi}
    
    and

    .. math::

        SMF = 1 - CMF
    
    `core_mass_fraction` is a mass fraction on the volatile-free, oxygen-free element basis; a true planetary CMF additionally requires the oxygen bound in mantle oxides and any light elements in the core. BSP and BP are normalized compositions, and so a uniform depletion of every element into the core (i.e., all alpha values are equal) is invisible to normalized compositions. This is the same convention used by Putirka & Rarick (2019), but it is an assumption, not something the compositions can test.

Mass retention factors
^^^^^^^^^^^^^^^^^^^^^^
    :math:`\alpha` and :math:`\phi` describe the partitioning of elements between mantle and core on a normalized concentration basis. They do not describe the mass of an element in either domain. For example, with :math:`\alpha_{Fe}` = 0.49 and an SMF = 0.9 (i.e., the mantle contains 90% of the planet's mass), only ~44% of the planet's Fe mass is in the mantle. Quoting alphas as retention fractions overstates mantle mass retention.

    If alpha and CMF or SMF are known, we can then calculate the mass fraction of any particular component in the mantle. We define a mass retention factor, phi, to describe the percentage of :math:`m_i` that resides in the planet's mantle.

    .. math::
        
        m_i^{mantle} = \alpha_{i} × SMF  

References
----------
Putirka, K.D. and Rarick, J.C. (2019) The Composition and Mineralogy of Rocky Exoplanets: A Survey of >4,000 Stars from the Hypatia Catalog