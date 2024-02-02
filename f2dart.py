#!/usr/bin/env python
from pathlib import Path
import re

def f2dart(glob: str):
    for path in Path('.').rglob(glob):
        print(path.resolve())
        lines = path.read_text().splitlines()
        start = True
        header_lines = []
        rest = False
        rest_lines = []
        for line in lines:
            if start and not line.strip():
                continue
            start = False

            if rest:
                rest_lines.append(line)
                continue

            if re.search('^[^\s]', line):
                header_lines.append(line)
            else:
                rest = True
                if line.strip():
                    rest_lines.append(line)

        with open(f'{str(path.resolve())[:-len(path.suffix)]}.header.txt', 'w') as header_file:
            header_file.write('\n'.join(header_lines))

        path.write_text('\n'.join(rest_lines))

f2dart('*.f')
f2dart('*.F')
f2dart('*.f90')
f2dart('*.F90')
