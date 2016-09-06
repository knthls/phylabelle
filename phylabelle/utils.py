# coding: utf-8
"""
little helpers
"""


class NoBooleanValueException(Exception):
    """
    The argument string doesn't hold a boolean value
    """
    pass


def parse_bool(b_str):
    """
    Parse the boolean value of a string. Parsing is case insensitive. Valid
    values are: ``1,0, true, false``
    :param str b_str:
    :raises NoBooleanValueException:
    :returns bool:
    """
    if b_str.lower() in ['1', 'true']:
        return True
    elif b_str.lower() in ['0', 'false']:
        return False
    else:
        raise NoBooleanValueException('The string \"{}\" doesn\'t '
                                      'seem to hold a boolean value'.format(b_str))
