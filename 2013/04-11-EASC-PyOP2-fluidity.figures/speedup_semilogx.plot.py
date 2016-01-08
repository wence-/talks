import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab

fig = pylab.figure('speedup_semilogx', figsize=(8, 6), dpi=300)
pylab.semilogx(np.load('elements.npy'), np.load('fluidity_speedup.npy'), 'k--o', lw=1, label='Fluidity (cores: 1)')
pylab.semilogx(np.load('elements.npy'), np.load('fluidity_mpi_speedup.npy'), 'g--s', lw=1, label='Fluidity (cores: 1)')
pylab.semilogx(np.load('elements.npy'), np.load('fluidity_pyop2_seq_speedup.npy'), 'r-d', lw=2, label='Fluidity-PyOP2 (backend: sequential)')
pylab.semilogx(np.load('elements.npy'), np.load('fluidity_pyop2_openmp_speedup.npy'), 'c-v', lw=2, label='Fluidity-PyOP2 (backend: openmp)')
pylab.semilogx(np.load('elements.npy'), np.load('fluidity_pyop2_mpi_speedup.npy'), 'm-<', lw=2, label='Fluidity-PyOP2 (backend: mpi)')
pylab.semilogx(np.load('elements.npy'), np.load('fluidity_pyop2_cuda_speedup.npy'), 'b-^', lw=2, label='Fluidity-PyOP2 (backend: cuda)')

pylab.legend(loc='lower right')
pylab.xlabel('Number of elements in the mesh')
pylab.ylabel('Relative speedup over Fluidity baseline')
pylab.title('Benchmark of an advection-diffusion problem for 100 time steps')
pylab.grid()
pylab.savefig('speedup_semilogx.svg', orientation='landscape', format='svg', transparent=True)
