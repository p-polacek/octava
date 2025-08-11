#!/usr/bin/env python


def var_csq(impact_tag=True):
    """Creates a list of tuples of variant consequences

    impact_tag: bool (default=True), add impact to display name

    returns: list, [(id_name, disp_name)]
    """
    csq_list = list()
    csq_dict = {"HIGH": "H", "MODERATE": "M", "LOW": "L", "MODIFIER": "-"}

    with open("static/variant_consequences.tsv", "r") as csq_in:
        for line in csq_in:
            id_name, disp_name, impact = line.strip().split("\t")[:3]
            if impact_tag:
                disp_name = f"{csq_dict[impact]} {disp_name}"
            csq_list.append((id_name, disp_name))

    return csq_list
