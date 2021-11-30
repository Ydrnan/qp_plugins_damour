#!/bin/python3

import sys

v1 = float(sys.argv[1])
v2 = float(sys.argv[2])

diff = abs(v1-v2)

if (diff > 10**(-10)):
    is_ok = False
else:
    is_ok = True

print("Diff: "+str(format(diff, ".3E")))
print("Valid: "+str(is_ok))

