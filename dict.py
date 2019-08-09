from typing import Any, Dict, List, Hashable
from iterable_wrappers import Quiet


def find_all_values(dicts: List[dict], key: Hashable, tqdm=Quiet) -> Dict[Any, int]:
    """
    Finds all values associated with a given key in a list of dictionaries.
    :param dicts: The list of dictionaries.
    :param key: The key to find values for.
    :param tqdm: The iterable wrapper to use.  Set to the tqdm class to have progress printed dynamically. Default
    Quiet, which does nothing.
    :return: A dictionary of (value=count) pairs.  Each value that is linked to the requested key in the dictionaries
    becomes a key of this return.  The associated value is how many times that key appeared as a value.  For instance,
    {"foo": 7, "bar": 12} indicates that the key had the value "foo" associated with it 7 times, and the value "bar"
    associated with 12 times.
    """
    ret = dict()

    for d in tqdm(dicts, desc="Sifting observations"):
        if key in d:
            val = d[key]
            try:
                ret[val] += 1
            except KeyError:
                ret[val] = 1

    return ret
