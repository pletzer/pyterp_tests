import numpy
from math import cos, sin, pi


def computeSphericalTriangleArea(lams, thes):
		
	eps = 1.e-10

	t0, t1, t2 = thes
	l0, l1, l2 = lams
	print 'thes = ', thes, ' lams = ', lams

	jac = (l1 - l0)*(t2 - t0) - (l2 - l0)*(t1 - t0)
	print 'jac = ', jac

	dt02 = t0 - t2
	dt01 = t0 - t1
	dt12 = t1 - t2
	area = 0.0

	cost0, cost1, cost2 = cos(t0), cos(t1), cos(t2)
	sint0, sint1 = sin(t0), sin(t1)

	if abs(dt02) < eps:
		# t0 == t2
		print 't0 == t2'
		area = jac*(dt01*sint0 + cost0 - cost1)/dt01**2
	else:
		if abs(dt12) < eps:
			# t1 == t2
			print 't1 == t2'
			area = jac*(-dt01*sint1 - cost0 + cost1)/dt01**2
			print 'area = ', area
		elif abs(dt01) < eps:
			print 't1 == t0'
			# t1 == t0
			area = jac*(dt02*sint0 + cost0 + cost2)/dt02**2
			print 'area = ', area
		else:
			print 'normal'
			# all the theta values are different
			dcos10 = (cost1 - cost0)/dt01
			dcos12 = (cost1 - cost2)/dt12
			area = jac*(dcos10 + dcos12)/dt02
	return area

#####################################################################

def test0():
	print computeSphericalTriangleArea([0., 2*pi, 2*pi], [0., 0., pi/2.])

def test1():
	print computeSphericalTriangleArea([0., 2*pi, pi], [0., 0., pi/2.])

def test2():
	area = computeSphericalTriangleArea([0., 2*pi, 0.], [0., 0., pi/2.]) + computeSphericalTriangleArea([2*pi, 2*pi, 0.], [0., pi/2., pi/2.])
	assert abs(area - 2*pi) < 1.e-10

def test3():
	area = 0.0
	area += computeSphericalTriangleArea([0., 2*pi, 0.], [0., 0., pi/2.]) + computeSphericalTriangleArea([2*pi, 2*pi, 0.], [0., pi/2., pi/2.])
	area -= computeSphericalTriangleArea([0., 2*pi, 0.], [0., 0., -pi/2.]) + computeSphericalTriangleArea([2*pi, 2*pi, 0.], [0., -pi/2., -pi/2.])
	assert abs(area - 4*pi) < 1.e-10	

def test4():
	area = 0.0
	area += computeSphericalTriangleArea([0., 2*pi, 0.], [0., 0., pi/2.]) + computeSphericalTriangleArea([2*pi, 2*pi, 0.], [0., pi/2., pi/2.])
	area += computeSphericalTriangleArea([0., 2*pi, 0.], [-pi/2., -pi/2., 0.0]) + computeSphericalTriangleArea([2*pi, 2*pi, 0.], [-pi/2., 0.0, 0.0])
	print 'total area should be 4*pi: ', area
	assert abs(area - 4*pi) < 1.e-10	

def test5():
	# half sphere
	# number of latitude and longitude sections
	nthe, nlam = 1, 8
	dthe, dlam = pi/float(2*nthe), 2*pi/float(nlam)
	totalArea = 0.0
	for i in range(nthe):
		the = 0.0 + i*dthe # top half
		for j in range(nlam):
			lam = 0.0 + j*dlam
			totalArea += computeSphericalTriangleArea([lam, lam + dlam, lam], [the, the, the + dthe])
			totalArea += computeSphericalTriangleArea([lam + dlam, lam + dlam, lam], [the, the + dthe, the + dthe])
	print 'total area = ', totalArea
	assert abs(totalArea - 2*pi) < 1.e-10

def test6():
	# half sphere
	# number of latitude and longitude sections
	nthe, nlam = 2, 1
	dthe, dlam = pi/float(2*nthe), 2*pi/float(nlam)
	totalArea = 0.0
	for i in range(nthe):
		the = 0.0 + i*dthe # top half
		for j in range(nlam):
			lam = 0.0 + j*dlam
			totalArea += computeSphericalTriangleArea([lam, lam + dlam, lam], [the, the, the + dthe])
			totalArea += computeSphericalTriangleArea([lam + dlam, lam + dlam, lam], [the, the + dthe, the + dthe])
	print 'total area = ', totalArea
	assert abs(totalArea - 2*pi) < 1.e-10



if __name__ == '__main__':
	"""
	test0()
	test1()
	test2()
	test3()
	test5()
	"""
	test6()
