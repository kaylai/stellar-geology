Definitions and Terms
=====================

This page is a work in progress.

The calculations implemented in `stellar_geology` are from Putirka & Rarick (2019) and benchmarked to the manuscript's supplementary data files.

Variable Definitions
--------------------
alphas
^^^^^^
    Alpha values as described in Putirka and Rarick (2019) represent the partitioning of elements between the bulk silicate portion of a planet and the core. This value can be set for any element known to `stellar_geology` For example:

    .. math::

        \alpha{Fe} = Fe^{BSP}/Fe^{BP}

    where :math:`Fe^{BSP}` is Fe in a bulk silicate planet (bulk planet, minus core) on a cation weight % basis (elemental weight proportions, absent anions) and :math:`Fe^{BP}` is the cation weight % of Fe for the bulk planet.

    Putirka & Rarick (2019) assume that :math:`\alpha_{Fe}^{Mercury}` = 0.0–0.12; :math:`\alpha_{Fe}^{Earth}` = 0.263–0.494; and :math:`\alpha_{Fe}^{Mars}` = 0.54–0.58.

References
----------
Putirka, K.D. and Rarick, J.C. (2019) The Composition and Mineralogy of Rocky Exoplanets: A Survey of >4,000 Stars from the Hypatia Catalog