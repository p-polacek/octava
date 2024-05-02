#!/usr/bin/env python
import json

"""
Haplosmasher v2.0 sample group manager
Author: Peter Poláček

Usage: To be called from within HaploSmasher v2.0 web interface
"""

def makeDeleteGroup(genome, data):
    """Creates or deletes a group of samples

    genome: str, name of genome saved in data/ directory
    data: dict(), info about the group to create or delete
    """
    json_file = f"data/{genome}/genome_data.json"

    with open(json_file, "r") as data_in:
        genome_data = json.load(data_in)
        genome_data["Groups"].setdefault(data["set"], list())

        # If making a group from a list of items
        if "name" in data.keys():
            group_entry = (data["items"], data["name"])

            # If a group with the given name exists, return error
            for group in genome_data["Groups"][data["set"]]:
                if group[1] == data["name"]:
                    return "Group with this name alrady exists."

            genome_data["Groups"][data["set"]].append(group_entry)

        # If deleting one or multiple groups
        else:
            filt = [group for group in genome_data["Groups"][data["set"]]
                    if group[1] not in data["items"]]
            genome_data["Groups"][data["set"]] = filt

    # Re-save the edited genome data
    with open(json_file, "w") as data_out:
        json.dump(genome_data, data_out)

    return None


if __name__ == "__main__":
    pass

