from xml.etree import ElementTree

def __tagged(name, value, attribute={}, indent=True):
    attribute_str = ''
    for key in attribute.keys():
        attribute_str += f' {key}="{attribute[key]}"'

    if value == '':
        return f'<{name}{attribute_str}/>'
    else:
        indented_value = ''
        for line in value.split('\n'):
            if line != '':
                indented_value += f'{"    " if indent else ""}{line}\n'
        header = f'<{name}{attribute_str}>\n'
        footer = f'</{name}>'
        return header + indented_value + footer

def read_xml(filename):
    tree = ElementTree.parse(filename)
    root = tree.getroot()
    name = root.tag
    return {
        name: __xml_to_dict(root)
    }

def __xml_to_dict(element):
    if len(element) == 0:
        if element.text is not None:
            element.attrib['text'] = element.text.strip()
        return element.attrib
    else:
        child_nodes = {child.tag: __xml_to_dict(child) for child in element}
        child_nodes.update(element.attrib)
        return child_nodes

def __dict_to_xml_helper(d, key_name, indent=True):
    attr_dict = {}
    text = ''
    for key in d.keys():
        child = d[key]
        if not isinstance(child, dict):
            if key == 'text':
                text = child
            else:
                attr_dict[key] = child
        else:
            text += __dict_to_xml_helper(child, key, indent=indent) + '\n'
    result = __tagged(key_name, text, attr_dict, indent=indent)
    return result

def dict_to_xml(d, indent=True):
    result = '<?xml version="1.0" encoding="UTF-8"?>\n'
    for key in d.keys():
        result += __dict_to_xml_helper(d[key], key, indent)
    return result
