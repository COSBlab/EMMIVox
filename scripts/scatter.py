import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
import sys

# read data from file
cc=[]; err=[]; dback=[]; dside=[]
for lines in open(sys.argv[1], "r").readlines():
    riga=lines.strip().split()
    # store stuff
    cc.append(float(riga[3]))
    # error
    err.append(float(riga[8]))
    # density backbone/sidechain
    dback.append(float(riga[9]))
    dside.append(float(riga[10]))

# do figure
fig, ax = plt.subplots()
ax.scatter(cc,err, s=dback, c=dside, alpha=0.9, cmap='coolwarm')

# labels
ax.set_xlabel('CC per residue', fontsize=15)
ax.set_ylabel('Relative error per residue', fontsize=15)

#ax.grid(True)
fig.tight_layout()

plt.show()
