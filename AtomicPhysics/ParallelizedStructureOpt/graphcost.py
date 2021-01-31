from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import autostructureoptimizer

## Do optimization
results = autostructureoptimizer.return_opt()


fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_title('6p and 6d Orbital Parameters versus. NIST Conf. Energy Difference')
ax.set_xlabel('6p Lambda')
ax.set_ylabel('6d Lambda')

ax.plot_surface(results[1], results[2], results[3])


ax.savefig('6p6dnistcompare.png')
