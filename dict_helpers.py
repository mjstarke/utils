from iterable_wrappers import Quiet
from operator import itemgetter
from typing import *


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


def pretty_print_dictionary(d: dict, print_percent: bool = True, print_total: bool = True,
                            total: Optional[float] = None, sorting: Union[None, str, Iterable] = "ka",
                            min_column_widths: Tuple[int, int, int] = (1, 1, 1),
                            compress_below: Optional[float] = None, column_separator: str = " | ") -> list:
    """
    Pretty-prints the contents of a dictionary.
    :param d: The dictionary to assess.
    :param print_percent: Whether to calculate percentages that each value makes up of the whole.  Fails if any values
    in d are not numeric.  Default True.
    :param print_total: Whether to print a row at the end for the sum of all the values.  Fails if any values in d are
    not numeric.  Default True.
    :param total: The number that will be printed if print_total or used to calculate percentage contributions if
    print_percent.  Default None, which automatically calculates total as the sum of dictionary values.
    :param sorting: If an iterable, specifies the order in which the keys will be printed.  If None, key ordering will
    be arbitrary (d.keys() is used).  If "ka", keys will be sorted ascending; "kd" sorts descending.  "va" will sort the
    keys by their corresponding values ascending; "vd" sorts by values descending.  Default "ka".
    :param min_column_widths: A triple of the minimum widths of the three printed columns (key, value, and percentage,
    respectively).  Default (1, 1, 1).
    :param compress_below: Items in d whose values are less than compress_below will be compressed into a single
    "(other)" entry at the bottom of the table (unless the sum of these items is 0).  Default None, which compresses no
    items.
    :param column_separator: The string which separates the printed columns.  Default ' | '.
    :return: None.  The results are printed in three columns: key, value, percentage (if do_percent).  The columns are
    automatically sized according to their contents.  If d contains no keys, nothing is printed.
    :raises TypeError: If sorting is a non-iterable non-string.
    :raises TypeError: If print_total or print_percent is true, but some values in the dictionary are not numeric.
    :raises ValueError: If min_column_widths does not have length 3.
    :raises ValueError: If any element of min_column_widths is less than 1.
    :raises ValueError: If sorting is an invalid string.
    :raises ZeroDivisionError: If total is zero (whether set explicitly or calculated automatically).
    :return: The list of keys in the order that they were printed, or an empty list if no keys were printed.  Note that
    the compression and total rows are not considered keys.
    """
    if len(min_column_widths) != 3:
        raise ValueError("Argument 'min_column_widths' must have length 3.")

    for width in min_column_widths:
        if width < 1:
            raise ValueError("All minimum column widths must be at least 1.")

    # If there are no keys in the dictionary, there's nothing to do.
    if len(d.keys()) == 0:
        return []

    # keys remains as this arbitrary list if sorting is None.
    keys = d.keys()

    if type(sorting) is str:
        if sorting == "ka":
            keys = sorted(d.keys())
        elif sorting == "kd":
            keys = sorted(d.keys())[::-1]
        elif sorting in "va":
            keys = [item[0] for item in sorted(d.items(), key=itemgetter(1))]
        elif sorting in "vd":
            keys = [item[0] for item in sorted(d.items(), key=itemgetter(1), reverse=True)]
        else:
            raise ValueError("Argument 'sorting' only accepts 'ka', 'kd', 'va', or 'vd' for strings.")
    elif sorting is not None:
        # If not None, sorting should be an iterable.  Test that here.
        try:
            _ = (e for e in d)
        except TypeError:
            raise TypeError("Argument 'sorting' must be None, 'ka', 'kd', 'va', 'vd', or an iterable.")
        keys = sorting

    # Use explicitly set total if available; otherwise, calculate total.
    if (total is None) and (print_percent or print_total):
        try:
            total = sum(d[key] for key in keys)
        except TypeError:
            raise TypeError("Total and percentages cannot be calculated with non-numeric values in the dictionary.")

    # If total is zero, we can't do percentage calculations.
    if print_percent and (total == 0):
        raise ZeroDivisionError("Total is 0; percentage contributions cannot be calculated.")

    # Chop off low values if requested.
    compressed = None
    if compress_below is not None:
        compressed = sum(d[key] for key in keys if d[key] < compress_below)
        keys = [key for key in keys if d[key] >= compress_below]

    # Generate contents of columns.
    column1 = [str(key) for key in keys]
    column2 = [str(d[key]) for key in keys]
    # Third column is blank if percent is not requested.
    column3 = ["{:7.2%}".format(d[key] / total) for key in keys] if print_percent else [""] * len(keys)

    if (compressed is not None) and (compressed > 0):
        column1.extend(["", "(other)"])
        column2.extend(["", str(compressed)])
        column3.extend(["", "{:7.2%}".format(compressed / total) if print_percent else ""])

    # Add two more rows for total if requested.
    if print_total:
        column1.extend(["", "(total)"])
        column2.extend(["", str(total)])
        column3.extend(["", "100.00%" if print_percent else ""])

    # Calculate the width that each column requires.
    column1_width = max(max(len(a) for a in column1), min_column_widths[0])
    column2_width = max(max(len(a) for a in column2), min_column_widths[1])
    column3_width = max(max(len(a) for a in column3), min_column_widths[2])

    # Create the format string with appropriate column width.
    fmt = "{:A}S{:>B}S{:C}".replace("A", str(column1_width)).replace("B", str(column2_width)).replace(
        "C", str(column3_width)).replace("S", column_separator)

    # Print each row.
    for a in range(len(column1)):
        print(fmt.format(column1[a], column2[a], column3[a]))

    return list(keys)
