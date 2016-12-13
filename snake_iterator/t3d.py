n0, n1, n2 = 3, 4, 5

def alternateDirection(iMinus1, n, i):
	return (1 - iMinus1 % 2)*i + (iMinus1 % 2)*(n - 1 - i)

inds = []
for i0 in range(n0):
	for i1 in range(n1):
		for i2 in range(n2):
			indexFlat = n1*n2*i0 + n2*i1 + i2
			#indexSnake = n1*j + (1 - j%2)*i + (n1 - 1 - i)*(j%2)
			indexSnake = n1*n2*i0 + \
			             alternateDirection(i0, n1, i1) * n2 + \
			             alternateDirection(i1, n2, i2) 

			inds.append(indexSnake)
			print('indexFlat = {} indexSnake = {}'.format(indexFlat, indexSnake))

inds.sort()
print(inds)

