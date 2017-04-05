from matplotlib import pylab

# tripolar to uniform 2560 x 5120 

nprocs = [1, 2, 4, 8, 16,]
times = [2.12e3, 1.18e3, 637, 368, 203]

pylab.plot(nprocs, [times[0]/t for t in times], 'ko', nprocs, [times[0]/t for t in times], 'b-')
pylab.title('ESMF conserve tripolar to uniform 2560x5120')
pylab.plot(nprocs, nprocs, 'k--')
pylab.xlabel('number of procs')
pylab.ylabel('speedup')
pylab.show()

