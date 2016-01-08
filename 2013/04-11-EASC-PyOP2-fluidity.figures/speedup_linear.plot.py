import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab

fig = pylab.figure('speedup_linear', figsize=(8, 6), dpi=300)
pylab.plot(np.load('elements.npy'), np.load('fluidity_speedup.npy'), 'k--o', lw=1, label='Fluidity+serial')
pylab.plot(np.load('elements.npy'), np.load('fluidity_mpi_speedup_fixed.npy'), 'g--s', lw=1, label='Fluidity+MPI (12 cores)')
pylab.plot(np.load('elements.npy'), np.load('fluidity_pyop2_seq_speedup.npy'), 'r-d', lw=2, label='PyOP2+serial')
# pylab.plot(np.load('elements.npy'), np.load('fluidity_pyop2_openmp_speedup.npy'), 'c-v', lw=2, label='Fluidity-PyOP2 (backend: openmp)')
pylab.plot(np.load('elements.npy'), np.load('fluidity_pyop2_cuda_speedup_fixed.npy'), 'b-^', lw=2, label='PyOP2+cuda')

pylab.legend(loc='right')
pylab.xlabel('Number of elements in the mesh')
pylab.ylabel('Relative speedup over Fluidity baseline')
pylab.title('Benchmark of an advection-diffusion problem for 100 time steps')
pylab.grid()
pylab.savefig('speedup_linear.pdf', orientation='landscape', format='pdf', transparent=True)
