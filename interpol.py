import numpy as np
import sys
import matplotlib.pyplot as plt
import math
import random as rnd
import argparse
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

verbose = 0

def linearInterpol(data):
    x = np.linspace(data[0][0],data[1][0],10000)
    y = []
    for val in x:
        y.append(data[0][1]+(data[1][1]-data[0][1])*(val-data[0][0])/(data[1][0]-data[0][0]))
    return x,y

def polyInterpol(data, degree, outline, past):
    coeff = degree+1
    rnd.seed(data[0][1])
    # Check if we can solve this using a Vandermonde matrix. 
    if(coeff < len(data)):
        cut = len(data)-coeff
        for i in range(0,cut):
            idx = rnd.randint(0,len(data)-1)
            data.remove(data[idx])
    matrix = []
    v = []
    for row in range(0,len(data)):
        matrix.append([])
        for col in range(0,coeff):
            matrix[row].append(math.pow(data[row][0], col))
        # Next prepare the goal vector
        v.append(data[row][1])
    # Solve matrix lambda = v
    lam = np.linalg.solve(matrix, v)
    x = np.linspace(0-past,data[len(data)-1][0],10000+outline)
    y = []
    for value in x:
        p = 0.0
        for i in range(len(lam)):
            p = p + lam[i]*math.pow(value,i)
        y.append(p)
    return x, y

# Use least squares to get a line
# See: https://neutrium.net/mathematics/least-squares-fitting-of-a-polynomial/
def LineOfBestFit(data, outline, past):
    x = [i[0] for i in data]
    y = [i[1] for i in data]
    meanX = 0.0
    for val in x:
        meanX = meanX+val
    meanX = meanX/len(x)
    meanY = 0.0
    for val in y:
        meanY = meanY+val
    meanY = meanY/len(y)
    m = 0.0
    up = 0.0
    down = 0.0
    for i in range(0,len(x)):
        up = up + (x[i]-meanX)*(y[i]-meanY)
        down = down+(x[i]-meanX)*(x[i]-meanX)
    m = up/down
    intercept = meanY-m*meanX
    
    dataX = np.linspace(data[0][0]-past,data[len(data)-1][0]+outline,10000)
    dataY = []
    for val in dataX:
        dataY.append(intercept+m*val)
    return dataX, dataY

# Use least squares to get a polynomial fit
# See: http://mathworld.wolfram.com/LeastSquaresFittingPolynomial.html
def PolyLeastSquares(data, degree, outline, past):
    x = [i[0] for i in data]
    y = [i[1] for i in data]
    matrix = []
    for row in range(0,degree+1):
        matrix.append([])
        for col in range(0,degree+1):
            value = 0.0
            for val in x:
                value = value + math.pow(val, col+row)
            matrix[row].append(value)
            
    v = []
    for row in range(0,degree+1):
        value = 0.0
        for i in range(0,len(x)):
            value = value + math.pow(x[i],row)*y[i]
        v.append(value)
    lam = np.linalg.solve(matrix, v)
    dataX = np.linspace(data[0][0]-past, data[len(data)-1][0]+outline,10000)
    dataY = []
    for value in dataX:
        p = 0.0
        for i in range(len(lam)):
            p = p + lam[i]*math.pow(value,i)
        dataY.append(p)
    if(verbose >  0):
        # Compute each f(x)
        newY = []
        for value in x:
            p = 0.0
            for i in range(len(lam)):
                p = p + lam[i]*math.pow(value,i)
            newY.append(p)
        # Get the sum and update the totalWeight
        currentWeight = 0.0
        for i in range(0,len(x)):
            value = newY[i]-y[i]
            value = value**2
            currentWeight = currentWeight + value
        print "Total error: ", currentWeight
        
    return dataX, dataY
    
# Use least squares to get an exponential fit
# See: http://mathworld.wolfram.com/LeastSquaresFittingExponential.html
def ExpLeastSquares(data, outline, past):
    x = [i[0] for i in data]
    y = [i[1] for i in data]
    matrix = []
    for row in range(0,2):
        matrix.append([])
        for col in range(0,2):
            value = 0.0
            for i in range(0,len(x)):
                if(row==0):
                    value = value + y[i]*math.pow(x[i],col)
                else:
                    value = value + y[i]*math.pow(x[i],col+1)
            matrix[row].append(value)
    v = []
    for row in range(0,2):
        value = 0.0
        for i in range(0,len(x)):
            value = value + y[i]*math.pow(x[i],row)*math.log(y[i])
        v.append(value)
    lam = np.linalg.solve(matrix, v)
    dataX = np.linspace(data[0][0]-past, data[len(data)-1][0]+outline,10000)
    dataY = []
    for value in dataX:
        p = 0.0
        p = math.exp(lam[0])*math.exp(lam[1]*value)
        dataY.append(p)
    return dataX, dataY

def mm(A, transposeA, B, transposeB):
    mA = []
    if(transposeA):
        for col in range(len(A[0])):
            mA.append([])
            for row in range(len(A)):
                mA[col].append(A[row][col])
    else:
        for row in range(len(A)):
            mA.append([])
            for col in range(len(A[row])):
                mA[row].append(A[row][col])
            
    mB = []
    if(transposeB):
        for col in range(len(B[0])):
            mB.append([])
            for row in range(len(B)):
                mB[col].append(B[row][col])
    else:
        for row in range(len(B)):
            mB.append([])
            for col in range(len(B[row])):
                mB[row].append(B[row][col])
    C = []
    for rowA in range(len(mA)):
        C.append([])
        for colB in range(len(mB[0])):
            value = 0.0
            for idx in range(len(mA[rowA])):
                value = value + mA[rowA][idx]*mB[idx][colB]
            C[rowA].append(value)
    return C

def mv(A, transposeA, V):
    mA = []
    if(transposeA):
        for col in range(len(A[0])):
            mA.append([])
            for row in range(len(A)):
                mA[col].append(A[row][col])
    else:
        for row in range(len(A)):
            mA.append([])
            for col in range(len(A[row])):
                mA[row].append(A[row][col])
    C = []
    for rowA in range(len(mA)):
        C.append([])
        value = 0.0
        for idx in range(len(mA[rowA])):
            value = value + mA[rowA][idx]*V[idx]
        C[rowA].append(value)
    return C

# Iteratively reqight least squares.
# Compute least squares fit.
# Compute weights
# Repeat until convergence
# https://cnx.org/contents/krkDdys0@12/Iterative-Reweighted-Least-Squ
def L1Norm(data, degree, outline, past):
    convergence = False
    x = [i[0] for i in data]
    y = [i[1] for i in data]
    weights = []
    rnd.seed(data[0][1])
    for i in range(0,len(x)):
        w =  rnd.random()
        weights.append(w)
    if(verbose > 1):
        print weights
    totalWeight = 999999999.999
    lam = []
    itera = 0
    while(not convergence):
        if(verbose > 0):
            print itera
            itera = itera + 1
        matrixX = []
        for row in range(0,len(x)):
            matrixX.append([])
            for col in range(0,degree+1):
                matrixX[row].append(math.pow(x[row], col))
        matrixXWeight = []
        for row in range(0,degree+1):
            matrixXWeight.append([])
            for col in range(0, len(x)):
                matrixXWeight[row].append(matrixX[col][row] * weights[col])
        # matrixXWeight is already transposed
        matrix = mm(matrixXWeight, False, matrixX, False)
        v = mv(matrixXWeight, False, y)
        lam = np.linalg.solve(matrix, v)
        # Update currentWeight and check if this is the same as before
        newY = []
        for value in x:
            p = 0.0
            for i in range(len(lam)):
                p = p + lam[i]*math.pow(value,i)
            newY.append(p)
        currentWeight = 0.0
        for i in range(0,len(x)):
            value = newY[i]-y[i]
            value = value**2
            currentWeight = currentWeight + value
        if(currentWeight == totalWeight):
            convergence = True
        if(currentWeight < totalWeight):
            totalWeight = currentWeight
        elif(verbose > 0):
            print "something fishy. currentWeight = ", currentWeight, " and totalWeight: ", totalWeight
        if(verbose > 1):
            print "currentWeight = ", currentWeight, " and totalWeight: ", totalWeight
        
    dataX = np.linspace(data[0][0]-past, data[len(data)-1][0]+outline,10000)
    dataY = []
    for value in dataX:
        p = 0.0
        for i in range(len(lam)):
            p = p + lam[i][0]*math.pow(value,i)
        dataY.append(p)
    return dataX, dataY

def MEstimators(data, outline, past):
    print "not yet implemented"

parser = ArgumentParser(
    description=
'''Runs several fits for different data. These fits are solutions
for the course 'Modelling I' at the Johannes Gutenberg-University in Mainz.''',
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-t', '--task', type=int, required=True, default=1,
                    help=
                    '''Choose the task (and therefore the data and fit) to solve.
Task 1: Linear interpolation with two datapoints.
Task 2: Polynomial interpolation by solving a Vandermonde matrix with degree 3 on 4 datapoints.
Task 3: Polynomial interpolation with least squares with degree 4 and 5 on 4 datapoints.
Task 4: Linear fitting with least squares on a dataset 'dataPCPC' or any
other specified with '-f'.
Task 5: Polynomial fitting with least squares on a dataset 'dataPCPC' or any
other specified with '-f'.
Task 6: Exponential fitting with least squares on a dataset 'dataPCPC' or any
other specified with '-f'.
Task 7: L1-Norm fitting on a dataset 'dataPCPC or any other specified with '-f'.
Here some datapoints are exchanged with outliers.
Task 8: M-Estimators approach to fit on a dataset 'dataPCPC or any other specified with '-f'.
Here some datapoints are exchanged with outliers.''')
parser.add_argument('-d', '--degree', type=int, required=False, default=3,
                    help='''Sets the degree to use for polynomial least square fitting''')
parser.add_argument('-f', '--datafile', type=str, required=False, default='dataPCPC',
                    help=
                    '''The datafile for tasks 4, 5 and 6. At default 'dataPCPC'
will be used. The file should have x-values at one line
and y-values on the next line. Then it can have x-values
again. The file won't be checked for errors.''')
parser.add_argument('-p', '--predict', type=int, required=False, default=0,
                    help=
                    '''The desired prediction for the dataset. Should be bigger than 0.''')
parser.add_argument('-past', type=int, required=False, default=0,
                    help=
                    '''The desired look into the past for the dataset. Should be bigger than 0.''')
parser.add_argument('-v', '--verbose', type=int, required=False, default=0,
                    help=
                    '''Set verbositiy.
0 = Minimal additional output.
1 = Some additional output like total error of fit.
2 = All additional output like current error of each iteration.''')
args = parser.parse_args()

degree = args.degree
task = args.task
datafile = args.datafile
outline = args.predict
verbose = args.verbose

if(outline < 0):
    print "Your desired prediction is smaller than 0. It is set to 0."
    outline = 0
past = args.past
if(past < 0):
    print "Your desired look into the past is smaller than 0. It is set to 0."
    past = 0

plotOne = True

if(task == 1):
    print "Task 1: Linear interpolation with few data"
    dataExpenses = [[7380,1420],[47640,21769]]
    xOriginal = [i[0] for i in dataExpenses]
    yOriginal = [i[1] for i in dataExpenses]
    x,y = linearInterpol(dataExpenses)
    figureName = "LinearInterpol"
    title = "Outcome depending on income"
    yTitle = "Outcome"
    xTitle = "Income (different countries)"
    grid = False
elif(task == 2):
    print "Task 2: Polynomial interpolation with degree 3"
    dataWeight = [[0,4.3],[10,39],[20,55],[30,80]]
    xOriginal = [i[0] for i in dataWeight]
    yOriginal = [i[1] for i in dataWeight]
    x,y = polyInterpol(dataWeight, 3, outline, past)
    figureName = "3Interpol"
    title = "Weight interpolation with degree 3"
    yTitle = "Weight in kg"
    xTitle = "Age in years"
    grid = False
elif(task == 3):
    print "Task 3: Polynomial interpolation with degree 4 and 5"
    dataWeight = [[0,4.3],[10,39],[20,55],[30,80]]
    xOriginal = [i[0] for i in dataWeight]
    yOriginal = [i[1] for i in dataWeight]
    x,y = PolyLeastSquares(dataWeight, 4, outline, past)
    x2Original = [i[0] for i in dataWeight]
    y2Original = [i[1] for i in dataWeight]
    x2,y2 = PolyLeastSquares(dataWeight, 5, outline, past)
    figureName = "45Interpol"
    title = "Weight interpolation with degree 4 and 5 (least squares)"
    yTitle = "Weight in kg"
    xTitle = "Age in years"
    grid = False
    plotOne = False
elif(task == 4):
    print "Task 4: Linear fitting and least squares"
    dataPCPC = []
    f = open(datafile, "r")
    year = True
    added = 0
    for line in f:
        if(year):
            year = False
            data = line.split()
            for x in data:
                dataPCPC.append([float(x)])
                added = added+1
        else:
            year = True
            data = line.split()
            for y in data:
                dataPCPC[len(dataPCPC)-added].append(float(y))
                added = added-1
    f.close()
    xOriginal = [i[0] for i in dataPCPC]
    yOriginal = [i[1] for i in dataPCPC]
    x,y = LineOfBestFit(dataPCPC, outline, past)
    figureName = "LinearLeastSquares"
    title = "PCPC in Germany with linear fitting and least squares"
    yTitle = "USD"
    xTitle = "Year"
    grid = False
elif(task == 5):
    print "Task 5: Polynomial fitting and least squares with degree ", degree
    dataPCPC = []
    f = open(datafile, "r")
    year = True
    added = 0
    for line in f:
        if(year):
            year = False
            data = line.split()
            for x in data:
                dataPCPC.append([float(x)])
                added = added+1
        else:
            year = True
            data = line.split()
            for y in data:
                dataPCPC[len(dataPCPC)-added].append(float(y))
                added = added-1
    f.close()
    xOriginal = [i[0] for i in dataPCPC]
    yOriginal = [i[1] for i in dataPCPC]
    if(degree >= len(xOriginal)):
        degree = len(xOriginal)-1
        print "Your given degree is too high for this few data. Degree is set to ", degree
    x,y = PolyLeastSquares(dataPCPC, degree, outline, past)
    figureName = "PolyLeastSquares"
    title = "PCPC in Germany (polynomial fitting, least squares, degree " + str(degree) + ")"
    yTitle = "USD"
    xTitle = "Year"
    grid = False
elif(task == 6):
    print "Task 6: Exponential fitting and least squares"
    dataPCPC = []
    f = open(datafile, "r")
    year = True
    added = 0
    for line in f:
        if(year):
            year = False
            data = line.split()
            for x in data:
                dataPCPC.append([float(x)])
                added = added+1
        else:
            year = True
            data = line.split()
            for y in data:
                dataPCPC[len(dataPCPC)-added].append(float(y))
                added = added-1
    f.close()
    xOriginal = [i[0] for i in dataPCPC]
    yOriginal = [i[1] for i in dataPCPC]
    x,y = ExpLeastSquares(dataPCPC, outline, past)
    figureName = "ExpLeastSquares"
    title = "PCPC in Germany with exponential fitting and least squares"
    yTitle = "USD"
    xTitle = "Year"
    grid = False
elif(task == 7):
    print "Task 7: L1-Norm Fitting with ", degree, " degrees"
    dataPCPC = []
    f = open(datafile, "r")
    year = True
    added = 0
    current = 0
    for line in f:
        if(year):
            year = False
            data = line.split()
            for x in data:
                dataPCPC.append([float(x)])
                added = added+1
        else:
            year = True
            data = line.split()
            for y in data:
                value = float(y)
                if(current == 1):
                    value = 1038
                if(current == 11):
                    value = 1369
                if(current == 29):
                    value = 192960
                if(current == 41):
                    value = 210120
                current = current + 1
                dataPCPC[len(dataPCPC)-added].append(float(y))
                added = added-1
    f.close()
    xOriginal = [i[0] for i in dataPCPC]
    yOriginal = [i[1] for i in dataPCPC]
    if(degree >= len(xOriginal)):
        degree = len(xOriginal)-1
        print "Your given degree is too high for this few data. Degree is set to ", degree
    x,y = L1Norm(dataPCPC, degree, outline, past)
    figureName = "L1Norm"
    title = "PCPC in Germany with outliers and L1-Norm fitting (" + degree + " degrees)"
    yTitle = "USD"
    xTitle = "Year"
    grid = False
elif(task == 8):
    print "Task 8: M-Estimators with ", degree, " degrees"
    dataPCPC = []
    f = open(datafile, "r")
    year = True
    added = 0
    current = 0
    for line in f:
        if(year):
            year = False
            data = line.split()
            for x in data:
                dataPCPC.append([float(x)])
                added = added+1
        else:
            year = True
            data = line.split()
            for y in data:
                value = float(y)
                if(current == 1):
                    value = 1038
                if(current == 11):
                    value = 1369
                if(current == 29):
                    value = 192960
                if(current == 41):
                    value = 210120
                current = current + 1
                dataPCPC[len(dataPCPC)-added].append(float(y))
                added = added-1
    f.close()
    xOriginal = [i[0] for i in dataPCPC]
    yOriginal = [i[1] for i in dataPCPC]
    x,y = MEstimators(dataPCPC, outline, past)
    figureName = "MEstimator"
    title = "PCPC in Germany with outliers and M-Estimators"
    yTitle = "USD"
    xTitle = "Year"
    grid = False
else:
    print "Error: No such task: ", task
    print "Use -h for help"
    sys.exit()

if(plotOne):
    plt.plot(x, y, xOriginal, yOriginal, 'ro')
    plt.xlabel(xTitle)
    plt.ylabel(yTitle)
    plt.title(title)
    plt.grid(grid)
    plt.savefig(figureName, dpi=600) 
else:
    plt.plot(x, y, 'r', x2, y2, 'b', xOriginal, yOriginal, 'ro', x2Original, y2Original, 'bo')
    plt.xlabel(xTitle)
    plt.ylabel(yTitle)
    plt.title(title)
    plt.grid(grid)
    plt.savefig(figureName, dpi=600) 