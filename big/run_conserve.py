from subprocess import call
import re
import argparse

parser = argparse.ArgumentParser(description='Exercise regridding')
parser.add_argument('--nprocs', type=int, dest='nprocs', default=1,
                    help='Number of procs (for programs supporting MPI execution)')
args = parser.parse_args()

def getEvaluationTime(filename):
	m = re.search(r'evaluation\s+([\d\.e\-]+)', open(filename, 'r').read())
	if m:
		return float(m.group(1))
	return None

def getWeightsTime(filename):
	m = re.search(r'weights\s+([\d\.e\-]+)', open(filename, 'r').read())
	if m:
		return float(m.group(1))
	return None

def getNumberOfInvalidPoints(filename):
	m = re.search(r'invalid points\:\s+(\d+)', open(filename, 'r').read())
	if m:
		return int(m.group(1))
	return None


dst_celldims = [(10, 20),
                (20, 40),
                (40, 80),
                (80, 160),
                (160, 320),
                (320, 640),
                (640, 1280),
                (1280, 2560),
                (2560, 5120),]
                #(5120, 10240)]


ns = []
esmf_eval = []
esmf_weights = []
esmf_eval_par = []
esmf_weights_par = []
for dstDims in dst_celldims:
	# generate the grids
	call(['python', 'generate_field.py', \
		'--dst_nj', '{}'.format(dstDims[0] + 1), \
		'--dst_ni', '{}'.format(dstDims[1] + 1), \
		])

	# comes from the coords_CF_ORCA12_GO6-2.nc file
        srcDims = (3606, 4322)

	srcN = srcDims[0] * srcDims[1]
	dstN = dstDims[0] * dstDims[1]
	ns.append(srcN * dstN)
	print('number of src * dst cells is {}'.format(srcN * dstN))

	# run esmf serial
	err = open('log.err', 'w')
	out = open('log.txt', 'w')
	call(['python', 'esmf_interp.py'], stdout=out, stderr=err)
	out.close()
	esmf_eval.append(getEvaluationTime('log.txt'))
	esmf_weights.append(getWeightsTime('log.txt'))

	if args.nprocs > 1:
	    # run esmf parallel
	    err = open('log.err', 'w')
	    out = open('log.txt', 'w')
	    call(['mpiexec', '-n', str(args.nprocs), 'python', 'esmf_conserve.py'], stdout=out, stderr=err)
            out.close()
	    esmf_eval_par.append(getEvaluationTime('log.txt'))
            esmf_weights_par.append(getWeightsTime('log.txt'))
            
        print('ns               = {}'.format(ns))
        print('esmf eval        = {}'.format(esmf_eval))
        print('esmf weights     = {}'.format(esmf_weights))
        print('esmf eval par    = {}'.format(esmf_eval_par))
        print('esmf weights par = {}'.format(esmf_weights_par))

# write to file
import re, time
ta = re.sub(' ', '_', time.asctime())
f = open('run_conserve-{}.csv'.format(ta), 'w')
if args.nprocs > 1:
    f.write('src_num_cells*dst_num_cells,esmf_eval,esmf_weights,esmf_eval_par,esmf_weights_par\n')
    for i in range(len(ns)):
	f.write('{},{},{},{},{}\n'.format(ns[i], esmf_eval[i], esmf_weights[i], esmf_eval_par[i], esmf_weights_par[i]))
else:
    f.write('src_num_cells*dst_num_cells,esmf_eval,esmf_weights\n')
    for i in range(len(ns)):
	f.write('{},{},{}\n'.format(ns[i], esmf_eval[i], esmf_weights[i]))
f.close()

from matplotlib import pylab
import matplotlib

legs = ['esmf eval', 'esmf wgts',]

pylab.loglog(ns, esmf_eval, 'ro', markersize=8)
pylab.loglog(ns, esmf_weights, 'rs', markersize=8)
if args.nprocs > 1:
    pylab.loglog(ns, esmf_eval_par, 'co', markersize=8) 
    pylab.loglog(ns, esmf_weights_par, 'cs', markersize=8)
    legs += ['esmf eval {}p'.format(args.nprocs), 'esmf wgts {}p'.format(args.nprocs)]
pylab.legend(legs, loc=4)

pylab.plot(ns, esmf_eval, 'r-')
pylab.plot(ns, esmf_weights, 'r--')
if args.nprocs > 1:
    pylab.plot(ns, esmf_eval_par, 'c-')
    pylab.plot(ns, esmf_weights_par, 'c--')

pylab.xlabel('num src cells * num dst cells')
pylab.ylabel('time [sec]')
pylab.title('conservative interpolation tripolar to rectilinear')
pylab.savefig('run_conserve.png')
pylab.show()
