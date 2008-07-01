#! /usr/bin/env python


IGNORE_KEYS = [
"Warning: Subprogram CSEND",
"Warning: Subprogram CRECV",
"Warning: Subprogram BCAST",
"Warning: Subprogram INDX1",
"Warning: Subprogram BYTE_WRITE",
"Warning: Subprogram DCLOCK",
"Warning: Subprogram Z_ROWSCALE",
"Warning: Subprogram Z_EXP",
"Warning: Subprogram MPI_",
"Warning: Subprogram DNCHBV",
"Warning: Subprogram ZBQLU01",
]



import sys

lines = sys.stdin.readlines()

i = 0

while i < len(lines):
    did_filtering = False
    for key in IGNORE_KEYS:
        if lines[i].startswith(key):
            i += 1
            while i < len(lines) and lines[i].startswith(" "):
                i += 1
            did_filtering = True
            break
    if not did_filtering:
        sys.stdout.write(lines[i])
        i += 1
