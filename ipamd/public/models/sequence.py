"""sequence class definition"""
from ipamd.public.utils.output import error

class Sequence:
    """base class for DNA, RNA and protein sequence"""
    unit = ''
    allowed_types = ''
    def __init__(self, name, sequence):
        self.__seq_name__ = name
        self.__sequence__ = self.__check(sequence)

    def __str__(self):
        return f"{self.__seq_name__}:{self.__sequence__}"

    def __len__(self):
        return len(self.__sequence__)

    def __getitem__(self, key):
        return self.__sequence__[key]

    def __iter__(self):
        return iter(self.__sequence__)

    def __check(self, sequence):
        if not isinstance(sequence, (str, list)):
            error("Sequence must be a string or list.")
            raise TypeError()
        result = ""
        for res in sequence:
            if res not in self.__class__.allowed_types:
                error(f"Invalid {self.__class__.unit} '{res}' in sequence.")
                raise ValueError()
            result += res
        return result

    @property
    def sequence(self): # pylint: disable=missing-function-docstring
        return self.__sequence__

    @property
    def name(self): # pylint: disable=missing-function-docstring
        return self.__seq_name__

    def __setitem__(self, key, value):
        value = self.__check(value)
        match key:
            case slice() | int():
                builder = list(self.__sequence__)
                builder[key] = value
                self.__sequence__ = "".join(builder)
            case str():
                self.__sequence__ = self.__sequence__.replace(key, value)
            case _:
                raise TypeError(f"Unsupported index type: {type(key)}")

    def __add__(self, other):
        match other:
            case Sequence():
                if self.__class__ != other.__class__:
                    error(f"Cannot add {self.__class__.__name__} and {other.__class__.__name__}.")
                    raise TypeError()
                new_sequence = self.__sequence__ + other.__sequence__
                return self.__class__(self.__seq_name__, new_sequence)
            case str():
                other = self.__check(other)
                return self.__class__(self.name, self.__sequence__ + other)
            case _:
                error(f"Unsupported type for addition: {type(other)}")
                raise TypeError()

class ProteinSequence(Sequence):
    """protein sequence class"""
    allowed_types = 'ACDEFGHIKLMNPQRSTVWYX'
    unit = 'amino acid'

class DNASequence(Sequence):
    """DNA sequence class"""
    allowed_types = 'ACGTX'
    unit = 'nucleotide'

class RNASequence(Sequence):
    """RNA sequence class"""
    allowed_types = 'ACGUX'
    unit = 'nucleotide'
