import re

def convert(bbox):
    if bbox[-1] in ['N', 'E']:
        multiplier = 1 
    else:
        multiplier = -1
    return multiplier * sum(float(x) / 60 ** n for n, x in enumerate(re.split('°|\'|\"', bbox[:-2])))