import numpy
from math import cos, sin, pi


def computeSphericalTriangleArea(lams, thes):
		
	eps = 1.e-10

	t0, t1, t2 = thes
	l0, l1, l2 = lams

	jac = (l1 - l0)*(t2 - t0) - (l2 - l0)*(t1 - t0)

	dt02 = t0 - t2
	dt01 = t0 - t1
	dt12 = t1 - t2
	area = 0.0

	cost0, cost1, cost2 = cos(t0), cos(t1), cos(t2)
	sint0, sint1 = sin(t0), sin(t1)

	if abs(dt02) < eps:
		area = jac*(dt01*sint0 + cost0 - cost1)/dt01**2
	else:
		if abs(dt12) < eps:
			# t1 == t2
			area = jac*(-dt01*sint1 - cost0 + cost1)/dt01**2
		elif abs(dt01) < eps:
			# t1 == t0
			area = jac*(dt02*sint0 + cost0 + cost2)/dt02**2
		else:
			# all the theta values are different
			dcos10 = (cost1 - cost0)/dt01
			dcos12 = (cost1 - cost2)/dt12
			area = jac*(dcos10 + dcos12)/dt02
	return area

#####################################################################

def test1():
	print computeSphericalTriangleArea([0., 2*pi, pi], [0., 0., pi/2.])

if __name__ == '__main__':
	test1()