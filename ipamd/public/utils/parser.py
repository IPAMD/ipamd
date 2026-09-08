from ipamd.public.utils.output import error
def value_of(value, lookup_table):
    if type(value) != str:
        return value
    else:
        if value.startswith('compute:'):
            substr = value[8:]
            index_of_open_curly = 0
            while index_of_open_curly != -1:
                index_of_open_curly = substr.find('{', index_of_open_curly)
                if index_of_open_curly == -1:
                    break
                else:
                    index_of_close_curly = substr.find('}', index_of_open_curly)
                    if index_of_close_curly == -1:
                        index_of_close_curly = len(substr)
                    variable_name = substr[index_of_open_curly + 1:index_of_close_curly]
                    try:
                        value = lookup_table[variable_name]
                    except KeyError:
                        error(f'Variable {variable_name} not found in lookup table.')
                        raise ValueError(f'Variable {variable_name} not found in lookup table.')
                    substr = substr[:index_of_open_curly] + str(value) + substr[index_of_close_curly + 1:]
            return eval(substr)
        elif value=='True':
            return True
        elif value=='False':
            return False
        else:
            try:
                return float(value)
            except ValueError:
                return value

def range_to_list(range_str):
    """
    convert range string to list
    """
    if range_str == '':
        return []
    ranges = range_str.split(',')
    result = []
    for r in ranges:
        if '-' in r:
            start, end = map(int, r.split('-'))
            result.extend(range(start, end+1))
        else:
            result.append(int(r))
    return result

def frame_index_list(target_frame):
    """
    convert target_frame (int, list, or range string) to a list of frame indices
    """
    if target_frame is None or target_frame == '':
        return []
    if isinstance(target_frame, str):
        return range_to_list(target_frame)
    if isinstance(target_frame, (list, tuple, range)):
        return list(target_frame)
    return [int(target_frame)]

def protein_range_split(protein_range):
    """
    split protein range into protein type and range
    """
    target_range = []
    target_type = ""
    if protein_range == '' or protein_range is None:
        pass
    elif '@' not in protein_range:
        target_type = protein_range
    else:
        target_type, target_range = protein_range.split('@')
        if target_range != 'all':
            target_range = range_to_list(target_range)
        else:
            target_range = []

    return target_type, target_range
