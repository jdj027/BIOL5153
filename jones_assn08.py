#! /usr/bin/env python

# Jake Jones
# Assignment 8

import re

# Name example
blurb = "Katherine went to the concert to see her favorite band, Catheryn and the the Cathryn's. She ran into her friend Kathryn, who introduced Katherine to her friend, Catherine. Together, they enjoyed the concert while texting inaudible snippets to their mutual friends, Kathrin and katharine."
# Debug list for testing
# names = "Katherine Catheryn Cathryn Kathryn Katherine Catherine Kathrin katharine"

# Regular expression to find all (K/C)ath(a/e/r/y/i/n/e)s
regularstring = '.atha?e?r?y?i?n?e?'

# Using re.search
print('Names found using re.search: ')
while True:
    match = re.search(regularstring, blurb, re.I)
    if not match:
        break
    print(match.group(0))
    blurb = blurb[match.end():]

# Redefine blurb
blurb = "Katherine went to the concert to see her favorite band, Catheryn and the the Cathryn's. She ran into her friend Kathryn, who introduced Katherine to her friend, Catherine. Together, they enjoyed the concert while texting inaudible snippets to their mutual friends, Kathrin and katharine."
print('\n')

# Using re.Findall instead
print('Names found using re.findall: ')
newmatch = re.findall(regularstring, blurb, re.I) # Findall waaaaay easier than re.search for this example
if newmatch:
    print(newmatch)
else:
    print('No matches.')
