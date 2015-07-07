#!/bin/env python

import sys
import os

print sys.argv[1:]

for fname in sys.argv[1:]:
    newname = fname.replace(" ", "_")
    newname = newname.replace("'", "")
    newname = newname.replace('"', "")
    newname = newname.replace('.fts', ".fits")
    
    os.rename(fname, newname)
