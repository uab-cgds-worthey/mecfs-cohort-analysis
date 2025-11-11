"""Compute distributions and spread of matches based on phenotype mixes of ME/CFS criteria to see
how that compares to phenotypic display of other diseases.

Per Carruthers et. al 2011 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3427890/):
A patient will meet the criteria for postexertional neuroimmune exhaustion (A), at least one
symptom from three neurological impairment categories (B), at least one symptom from three
immune/gastro-intestinal/genitourinary impairment categories (C), and at least one symptom from
energy metabolism/transport impairments (D).

The goal of this script is to compute the possible mix of phenotype terms allowing assignments
from the above category and highlighting how dissimilar patients meeting these criteria could be.
"""

import collections
import itertools
import six


def load_phenos_by_group(pheno_file):
    """Convert TSV file with disease phenotypes into a dictionary of phenotypes per
    group in the third column of the file.

    Args:
        pheno_file (str): File path to the disease phenotypes TSV file

    Returns:
        dict: key = third column group, value = list of phenotype terms for that category
    """
    group_dict = dict()
    with open(pheno_file, "rt") as phenofp:
        for ln_index, line in enumerate(phenofp):
            # skip first line as that's just the header
            if ln_index == 0:
                continue
            cols = line.strip().split("\t")
            if cols[2] not in group_dict:
                group_dict[cols[2]] = []

            group_dict[cols[2]].append(cols[0])

    return group_dict


def iterable(arg):
    # https://stackoverflow.com/a/44328500/2892199
    return isinstance(arg, collections.abc.Iterable) and not isinstance(
        arg, six.string_types
    )


def get_mecfs_pheno_combos(pheno_groups_dict):
    """Create all possible combinations of MECFS phenotypes that minimally satisfy
    the Carruthers et. al 2011 MECFS diagnostic criteria. The result will not
    include PENE phenotypes as those are manditory diagnostic criteria for all
    MECFS cases and can be added after the fact.

    From Carruthers et. al 2011:
    A patient will meet the criteria for postexertional neuroimmune exhaustion (A), at least one
    symptom from three neurological impairment categories (B), at least one symptom from three
    immune/gastro-intestinal/genitourinary impairment categories (C), and at least one symptom
    from energy metabolism/transport impairments (D).

    Args:
        pheno_groups_dict (_type_): Dict with keys corresponding to the symptom subsections
        outlined in Carruthers et. al 2011

    Returns:
        List[List[str]]: list of MECFS minimal phenotypic avatars (as a list of HPO IDs)
    """
    # From Carruthers et. al 2011:
    # A patient will meet the criteria for postexertional neuroimmune exhaustion (A), at least one
    # symptom from three neurological impairment categories (B), at least one symptom from three
    # immune/gastro-intestinal/genitourinary impairment categories (C), and at least one symptom
    # from energy metabolism/transport impairments (D).

    # compute all combinations of phenotypes that satisfy category B
    group_b_phenos = [
        pheno_groups_dict["Neurocognitive"],
        pheno_groups_dict["Pain"],
        pheno_groups_dict["Sleep"],
        pheno_groups_dict["Neurosensory"],
    ]
    # create list of unique combinations of the groups with category B (mixes of the 3 out of 4 groups possible)
    group_combos = list(itertools.combinations(group_b_phenos, 3))
    # create all unique combinations a single phenotype term from group combos (i.e., unique phenotype
    # combinations of 3 terms where each term comes from 1 of the 3 groups)
    cat_b_combos = []
    for groupcombo in group_combos:
        cat_b_combos.extend(itertools.product(*groupcombo))

    # compute all combinations of phenotypes that satisfy category C
    cat_c_combos = list(
        itertools.combinations(pheno_groups_dict["Immune-Intestinal"], 3)
    )

    # combine all phenotype categories to build all possible phenotypic profiles
    pheno_product = itertools.product(
        pheno_groups_dict["Energy"],
        cat_b_combos,
        cat_c_combos,
    )
    # flatten the tuples produced for category B and C phenotypes so that profiles are lists of phenotype strings
    phenotypic_profiles = []
    for ppro in pheno_product:
        tmplist = []
        for prow in ppro:
            if iterable(prow):
                tmplist.extend(prow)
            else:
                tmplist.append(prow)
        phenotypic_profiles.append(tmplist)

    return phenotypic_profiles
