"""Common helpers for checking if command line arguments are valid"""

import os


def is_valid_file(parser, arg):
    """checks that a file path is valid

    Parameters
    ----------
    parser : ArgParse.ArgumentParser
        the argument parser object
    arg : str
        file path to check

    Returns
    -------
    arg : str
        the path to the file if valid
    """
    if not os.path.isfile(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg  # return the arg


def is_positive_integer(parser, arg):
    """checks that an argument is a positive int

    Parameters
    ----------
    parser : ArgParse.ArgumentParser
        the argument parser object
    arg : str
        value to check

    Returns
    -------
    ivalue : int
        integer representation of `arg`
    """
    ivalue = int(arg)
    if ivalue < 0:
        raise parser.ArgumentTypeError("%s is an invalid positive int value" % arg)
    else:
        return ivalue


def is_valid_path(parser, arg):
    """checks that folder path is valid

    Parameters
    ----------
    parser : ArgParse.ArgumentParser
        the argument parser object
    arg : str
        file path to check

    Returns
    -------
    arg : str
        the str path for the folder if valid

    """
    if not os.path.isdir(arg):
        parser.error("The path %s does not exist!" % arg)
    else:
        return arg  # return the arg


