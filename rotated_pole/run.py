from subprocess import call
import re
import argparse

parser = argparse.ArgumentParser(description='Exercise regridding')
parser.add_argument('--nprocs', type=int, dest='nprocs', default=1,
                    help='Number of procs (for programs supporting MPI execution)')

args = parser.parse_args()

def getEvaluationTime(filename):
	m = re.search(r'evaluation\s+([\d\.e\-\+]+)', open(filename, 'r').read())
	if m:
		return float(m.group(1))
	return None

def getWeightsTime(filename):
	m = re.search(r'weights\s+([\d\.e\-\+]+)', open(filename, 'r').read())
	if m:
		return float(m.group(1))
	return None

def getNumberOfInvalidPoints(filename):
	m = re.search(r'invalid points\:\s+(\d+)', open(filename, 'r').read())
	if m:
		return int(m.group(1))
	return None


src_celldims = [(10, 20),
                (20, 40),
                (40, 80),]
                #(80, 160),
                #(160, 320),
                #(320, 640),
                #(640, 1280),
                #(1280, 2560),]
                #(2560, 5120),
                #(5120, 10240)]


ns = []
libcf_interp_eval = []
libcf_interp_weights = []
esmf_interp_eval = []
esmf_interp_weights = []
esmf_conserve_eval = []
esmf_conserve_weights = []
esmf_interp_eval_par = []
esmf_interp_weights_par = []
esmf_conserve_eval_par = []
esmf_conserve_weights_par = []

for srcDims in src_celldims:
	# generate the grids
	dstDims = (srcDims[0]//2, srcDims[1]//2)
	call(['python', 'generate_field.py', \
		'--src_nj', '{}'.format(srcDims[0] + 1), \
		'--src_ni', '{}'.format(srcDims[1] + 1), \
		'--dst_nj', '{}'.format(dstDims[0] + 1), \
		'--dst_ni', '{}'.format(dstDims[1] + 1), \
		])

	srcN = srcDims[0] * srcDims[1]
	dstN = dstDims[0] * dstDims[1]
	ns.append(srcN * dstN)
	print('number of src * dst cells is {}'.format(srcN * dstN))

	# run esmf bilinear (serial)
	err = open('log.err', 'w')
	out = open('log.txt', 'w')
	call(['python', 'esmf_interp.py'], stdout=out, stderr=err)
	out.close()
	esmf_interp_eval.append(getEvaluationTime('log.txt'))
	esmf_interp_weights.append(getWeightsTime('log.txt'))

	if args.nprocs > 1:
		# run esmf conserve (parallel)
		err = open('log.err', 'w')
		out = open('log.txt', 'w')
		call(['mpiexec', '-n', str(args.nprocs), 'python', 'esmf_conserve.py'], stdout=out, stderr=err)
		out.close()
		esmf_conserve_eval_par.append(getEvaluationTime('log.txt'))
		esmf_conserve_weights_par.append(getWeightsTime('log.txt'))

		# run esmf bilinear (parallel)
		err = open('log.err', 'w')
		out = open('log.txt', 'w')
		call(['mpiexec', '-n', str(args.nprocs), 'python', 'esmf_interp.py'], stdout=out, stderr=err)
		out.close()
		esmf_interp_eval_par.append(getEvaluationTime('log.txt'))
		esmf_interp_weights_par.append(getWeightsTime('log.txt'))

	# run esmf conserve (serial)
	err = open('log.err', 'w')
	out = open('log.txt', 'w')
	call(['python', 'esmf_conserve.py'], stdout=out, stderr=err)
	out.close()
	esmf_conserve_eval.append(getEvaluationTime('log.txt'))
	esmf_conserve_weights.append(getWeightsTime('log.txt'))

	# run libcf (bilinear)
	err = open('log.err', 'w')
	out = open('log.txt', 'w')
	call(['python', 'libcf_interp.py'], stdout=out, stderr=err)
	out.close()
	libcf_interp_eval.append(getEvaluationTime('log.txt'))
	libcf_interp_weights.append(getWeightsTime('log.txt'))
	numFails = getNumberOfInvalidPoints('log.txt')
	if numFails != 0:
		print('*** {} libcf interp failures'.format(numFails))

	print(ns)
	print(esmf_interp_eval)
	print(esmf_interp_weights)
	print(esmf_conserve_eval)
	print(esmf_conserve_weights)
	print(libcf_interp_eval)
	print(libcf_interp_weights)

# write to file
import re, time
ta = re.sub(' ', '_', time.asctime())
f = open('run_node_interp-{}.csv'.format(ta), 'w')
f.write('src_num_cells*dst_num_cells,esmf_interp_eval,esmf_interp_weights,esmf_conserve_eval,esmf_conserve_weights,libcf_interp_eval,libcf_interp_weights\n')
for i in range(len(ns)):
	f.write('{},{},{},{},{},{},{}\n'.format(ns[i], esmf_interp_eval[i], esmf_interp_weights[i], libcf_interp_eval[i], libcf_interp_weights[i], esmf_conserve_eval[i], esmf_conserve_weights[i]))
f.close()

from matplotlib import pylab
import matplotlib

legs = []

legs.append('esmf lin eval')
pylab.loglog(ns, esmf_interp_eval, 'r-')
legs.append('esmf lin wgts')
pylab.loglog(ns, esmf_interp_weights, 'r--')
legs.append('libcf lin eval')
pylab.loglog(ns, libcf_interp_eval, 'b-')
legs.append('libcf lin wgts')
pylab.loglog(ns, libcf_interp_weights, 'b--')
legs.append('esmf con eval')
pylab.loglog(ns, esmf_conserve_eval, 'm-')
legs.append('esmf con wgts')
pylab.loglog(ns, esmf_conserve_weights, 'm--')
if args.nprocs > 1:
	legs.append('esmf lin eval {}p'.format(args.nprocs))
	pylab.loglog(ns, esmf_interp_eval_par, 'r-', linewidth=2)
	legs.append('esmf lin wgts {}p'.format(args.nprocs))
	pylab.loglog(ns, esmf_interp_weights_par, 'r--', linewidth=2)
	legs.append('esmf con eval {}p'.format(args.nprocs))
	pylab.loglog(ns, esmf_conserve_eval_par, 'm-', linewidth=2)
	legs.append('esmf con wgts {}p'.format(args.nprocs))
	pylab.loglog(ns, esmf_conserve_weights_par, 'm--', linewidth=2)

pylab.legend(legs, loc=2)
pylab.xlabel('num src cells * num dst cells')
pylab.ylabel('time [sec]')
pylab.title('regridding rotated pole to rectilinear')
pylab.savefig('run.png')
#pylab.show()