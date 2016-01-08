import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab

fig = pylab.figure('runtime_semilogx', figsize=(8, 6), dpi=300)
pylab.semilogx(np.load('elements.npy'), np.load('fluidity.npy'), 'k--o', lw=1, label='Fluidity (cores: 1)')
pylab.semilogx(np.load('elements.npy'), np.load('fluidity_mpi.npy'), 'g--s', lw=1, label='Fluidity (cores: 1)')
pylab.semilogx(np.load('elements.npy'), np.load('fluidity_pyop2_seq.npy'), 'r-d', lw=2, label='Fluidity-PyOP2 (backend: sequential)')
pylab.semilogx(np.load('elements.npy'), np.load('fluidity_pyop2_openmp.npy'), 'c-v', lw=2, label='Fluidity-PyOP2 (backend: openmp)')
pylab.semilogx(np.load('elements.npy'), np.load('fluidity_pyop2_mpi.npy'), 'm-<', lw=2, label='Fluidity-PyOP2 (backend: mpi)')
pylab.semilogx(np.load('elements.npy'), np.load('fluidity_pyop2_cuda.npy'), 'b-^', lw=2, label='Fluidity-PyOP2 (backend: cuda)')

pylab.legend(loc='upper left')
pylab.xlabel('Number of elements in the mesh')
pylab.ylabel('Overall runtime in seconds')
pylab.title('Benchmark of an advection-diffusion problem for 100 time steps')
pylab.grid()
pylab.savefig('runtime_semilogx.svg', orientation='landscape', format='svg', transparent=True)
