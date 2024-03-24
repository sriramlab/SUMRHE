import numpy as np
def partition_string(s, part):
    size = len(s)//part
    return [s[i:i+size] for i in range(0, len(s), size)]


l = np.arange(142)
p = partition_string(l, 10)
print(p)
print(len(p))