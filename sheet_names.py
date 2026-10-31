"""Excel sheet names for the comparisons, shared by the DESeq2 and TF-activity workbooks so their
sheets match. Excel allows 31 characters per sheet name."""
import re


def sheet_name(d, c):
    return f"{d['sheet']} " + re.sub(r", all \d+ pairs", "", c["label"]).replace(", hormone-naive", " naive")


def check(names):
    assert all(len(n) <= 31 for n in names) and len(set(names)) == len(names), names
