from subprocess import call
import re

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
                (40, 80),
                (80, 160),
                (160, 320),
                (320, 640),
                (640, 1280),
                (1280, 2560),]
                #(2560, 5120),
                #(5120, 10240)]


ns = []
libcf_interp_eval = []
libcf_interp_interp_weights = []
esmf_interp_eval = []
esmf_interp_weights = []
esmf_conserve_eval = []
esmf_conserve_weights = []

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

	# run esmf bilinear
	err = open('log.err', 'w')
	out = open('log.txt', 'w')
	call(['python', 'esmf_interp.py'], stdout=out, stderr=err)
	out.close()
	esmf_interp_eval.append(getEvaluationTime('log.txt'))
	esmf_interp_weights.append(getWeightsTime('log.txt'))

	# run esmf conserve
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
	f.write('{},{},{},{},{}\n'.format(ns[i], esmf_interp_eval[i], esmf_interp_interp_interp_weights[i], libcf_interp_interp_eval[i], libcf_interp_weights[i]))
f.close()

from matplotlib import pylab
import matplotlib
pylab.loglog(ns, esmf_interp_interp_interp_interp_eval, 'ro', markersize=8)
pylab.loglog(ns, esmf_interp_interp_interp_weights, 'rs', markersize=8)
pylab.loglog(ns, libcf_interp_interp_eval, 'bo', markersize=8) 
pylab.loglog(ns, libcf_interp_weights, 'bs', markersize=8)
pylab.legend(['esmf lin eval', 'esmf lin wgts', 'libcf lin eval', 'libcf lin wgts'], loc=2)
pylab.plot(ns, libcf_interp_interp_interp_interp_eval, 'b--', ns, libcf_interp_interp_interp_weights, 'b--', \
	       ns, esmf_interp_interp_eval, 'r--', ns, esmf_interp_weights, 'r--')
pylab.xlabel('num src cells * num dst cells')
pylab.ylabel('time [sec]')
pylab.title('bilinear interpolation')
pylab.show()