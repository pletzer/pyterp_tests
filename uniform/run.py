from subprocess import call
import re
import numpy

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

def getInterpTime(filename):
        m = re.search(r'interp\s+([\d\.e\-]+)', open(filename, 'r').read())
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

sigrid_conserve_eval = []
sigrid_conserve_weights = []

iris_interp = []
iris_conserve = []

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

	# run esmf linear
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

	# run libcf
	err = open('log.err', 'w')
	out = open('log.txt', 'w')
	call(['python', 'libcf_interp.py'], stdout=out, stderr=err)
	out.close()
	libcf_interp_eval.append(getEvaluationTime('log.txt'))
	libcf_interp_weights.append(getWeightsTime('log.txt'))

	# run sigrid
	err = open('log.err', 'w')
	out = open('log.txt', 'w')
	call(['python', 'sigrid_conserve.py'], stdout=out, stderr=err)
	out.close()
	sigrid_conserve_eval.append(getEvaluationTime('log.txt'))
	sigrid_conserve_weights.append(getWeightsTime('log.txt'))

	# run iris linear
	err = open('log.err', 'w')
	out = open('log.txt', 'w')
	call(['python', 'iris_interp.py'], stdout=out, stderr=err)
	out.close()
	iris_interp.append(getInterpTime('log.txt'))

	# run iris conserve
	err = open('log.err', 'w')
	out = open('log.txt', 'w')
	call(['python', 'iris_conserve.py'], stdout=out, stderr=err)
	out.close()
	iris_conserve.append(getInterpTime('log.txt'))

print('ns = {}'.format(ns))
print('esmf_interp_eval = {}'.format(esmf_interp_eval))
print('emsf_interp_weights = {}'.format(esmf_interp_weights))
print('libcf_interp_eval = {}'.format(libcf_interp_eval))
print('libcf_interp_weights = {}'.format(libcf_interp_weights))
print('iris_interp = {}'.format(iris_interp))
print('emsf_conserve_eval = {}'.format(esmf_conserve_eval))
print('esmf_conserve_weights = {}'.format(esmf_conserve_weights))
print('sigrid_conserve_eval = {}'.format(sigrid_conserve_eval))
print('sigrid_conserve_weights = {}'.format(sigrid_conserve_weights))
print('iris_conserve = {}'.format(iris_conserve))

# write to file
import re, time
ta = re.sub(' ', '_', time.asctime())
f = open('run_node_interp-{}.csv'.format(ta), 'w')
f.write(\
'src_num_cells*dst_num_cells,esmf_interp_eval,esmf_interp_weights,libcf_interp_eval,libcf_interp_weights,iris_interp,emsf_conserve_eval,esmf_conserve_weights,sigrid_conserve_eval,sigrid_conserve_weights,iris_conserve\n')
for i in range(len(ns)):
	f.write('{},{},{},{},{},{},{},{},{},{},{}\n'.format(ns[i], esmf_interp_eval[i], esmf_interp_weights[i], libcf_interp_eval[i], libcf_interp_weights[i], iris_interp[i], iris_interp[i],\
                                                 esmf_conserve_eval[i], esmf_conserve_weights[i], sigrid_conserve_eval[i], sigrid_conserve_weights[i], iris_conserve[i]))
f.close()

from matplotlib import pylab
import matplotlib
pylab.loglog(ns, numpy.array(esmf_interp_eval) + numpy.array(esmf_interp_weights), 'r--')
pylab.loglog(ns, numpy.array(libcf_interp_eval) + numpy.array(libcf_interp_weights), 'b--')
pylab.loglog(ns, iris_interp, 'g--')
pylab.loglog(ns, numpy.array(esmf_conserve_eval) + numpy.array(esmf_conserve_weights), 'r-')
pylab.loglog(ns, numpy.array(sigrid_conserve_eval) + numpy.array(sigrid_conserve_weights), 'c-')
pylab.loglog(ns, iris_conserve, 'g-')
legs = ['esmf lin', 'libcf lin', 'iris lin', 'esmf con', 'sigrid con', 'iris con']
pylab.legend(legs, loc=2)

pylab.xlabel('num src cells * num dst cells')
pylab.ylabel('time [sec]')
pylab.title('regridding (uni -> uni)')
pylab.show()
