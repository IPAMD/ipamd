# pylint: disable=missing-module-docstring
# pylint: disable=missing-function-docstring
import rich
from ipamd.public.models.data import String

def func(ref: String, target: String, print_ref=True, seperate_with_bracket=True):
    """
    Print the difference between two sequences
    :param seq1: the first sequence
    :param seq2: the second sequence
    :param print_ref: whether to print the reference sequence
    """
    type_ = ""
    ref_data = ref.data
    target_data = target.data
    for char1, char2 in zip(ref_data, target_data):
        if char1 == char2:
            type_ += " "
        elif char1 == "-" or char2 == "-":
            type_ += "d"
        else:
            type_ += "r"

    format_ = []
    format_entry = {}
    last_type = ""
    for i, char in enumerate(type_):
        if char != last_type:
            if last_type == "r" and char == "d" and seperate_with_bracket:
                continue
            if format_entry:
                format_entry["end"] = i - 1
                format_.append(format_entry)
                format_entry = {}
            format_entry["start"] = i
            format_entry["type"] = char
            last_type = char
    format_entry["end"] = len(type_) - 1
    format_.append(format_entry)

    def print_seq_diff(seq, format_):
        for format_entry in format_:
            type_ = format_entry["type"]
            start = format_entry["start"]
            end = format_entry["end"]
            if type_ == "r":
                template = "[{}]" if seperate_with_bracket else "{}"
                rich.print(
                    f"[yellow]{template.format(seq[start: end + 1])}[/yellow]",
                    end=""
                )
            elif type_ == "d":
                rich.print(f"[red]{seq[start: end + 1]}[/red]", end="")
            else:
                rich.print(f"[green]{seq[start: end + 1]}[/green]", end="")
        rich.print()

    if print_ref:
        print_seq_diff(ref_data, format_)
    print_seq_diff(target_data, format_)
