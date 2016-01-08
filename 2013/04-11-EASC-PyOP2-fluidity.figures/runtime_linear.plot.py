import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab

fig = pylab.figure('runtime_linear', figsize=(8, 6), dpi=300)
pylab.plot(np.load('elements.npy'), np.load('fluidity.npy'), 'k--', lw=2, label='Fluidity (cores: 1)')
pylab.plot(np.load('elements.npy'), np.load('fluidity_mpi_fixed.npy'), 'g--', lw=2, label='Fluidity (cores: 12)')
pylab.plot(np.load('elements.npy'), np.load('fluidity_pyop2_seq.npy'), 'r-', lw=2, label='PyOP2+serial')
pylab.plot(np.load('elements.npy'), np.load('fluidity_pyop2_openmp.npy'), 'c-', lw=2, label='PyOP2+openmp')
pylab.plot(np.load('elements.npy'), np.load('fluidity_pyop2_mpi.npy'), 'm-', lw=2, label='PyOP2+mpi')
pylab.plot(np.load('elements.npy'), np.load('fluidity_pyop2_cuda.npy'), 'b-', lw=2, label='PyOP2+cuda')

pylab.legend(loc='upper left')
pylab.xlabel('Number of elements in the mesh')
pylab.ylabel('Overall runtime in seconds')
pylab.grid()
pylab.savefig('runtime_linear.pdf', orientation='landscape', format='pdf', transparent=True)
