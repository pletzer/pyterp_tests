n0, n1 = 3, 4

inds = []
for j in range(n0):
	for i in range(n1):
		indexFlat = n1*j + i
		indexSnake = n1*j + (1 - j%2)*i + (n1 - 1 - i)*(j%2)
		inds.append(indexSnake)
		print('indexFlat = {} indexSnake = {}'.format(indexFlat, indexSnake))

inds.sort()
print(inds)

