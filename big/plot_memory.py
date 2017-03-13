import pandas
from matplotlib import pylab

data = pandas.read_csv('memory_results.csv')
totNumCells = 3606 * 4322 * data['Dst grid points']
memPeak = data[' MaxRSS (kbytes)'] / 1.e6
memTop = data[' top memory GB']

print memPeak

pylab.plot(totNumCells, memPeak, 'm-' , totNumCells, memTop, 'c-')
pylab.legend(['peak from SLURM', 'peak from top'], loc=2)
pylab.plot(totNumCells, memPeak, 'ks' , totNumCells, memTop, 'ks')
pylab.xlabel('num src cells * num dst cells')
pylab.ylabel('Gigabytes')
pylab.title('Memory consumption')
pylab.show()