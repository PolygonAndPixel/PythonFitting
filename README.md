# PythonFitting
In this script I am using different calculations to fit some data to a function.
These fits are solutions some tasks of the course **Modellierung I** by Prof. Michael Wand
at the Johannes Gutenberg-University in Mainz. 
I give no guarantee for correctness, style or efficiency (which is obvious if 
you look at my matrix-matrix multiplication for instance). However if you want
to see some basic implementations, feel free to take a look, copy some code and
play with it.
If you want to get to know something about style, please visit 
[PEP8](https://www.python.org/dev/peps/pep-0008/) or use at least
[pylint](https://www.pylint.org/) which has some nice GUI interfaces like 
[Spyder IDE](https://github.com/spyder-ide/spyder).
Do not hesitate to contact me if you find any mistakes or simply open an issue 
or make a pull request.

## Usage
Simply type

`python interpol.py [-t task]`

to get started. There are different options for most plots:

| Options | Description |
|---------|-------------|
| '-t', '--task' | type=int, required=True, default=1, choose the task (and therefore the data and fit) to solve.|
|'-d', '--degree'| type=int, required=False, default=3, sets the degree to use for polynomial least square fitting.|
|'-f', '--datafile'| type=str, required=False, default='dataPCPC', the datafile for tasks 4, 5 and 6. |
|'-p', '--predict'| type=int, required=False, default=0, the desired prediction for the dataset. Should be bigger than 0.|
|'-past', type=int| required=False, default=0, the desired look into the past for the dataset. Should be bigger than 0.|
|'-v', '--verbose'| type=int, required=False, default=0, set verbositiy from 0 to 2.|
|'-h', '--help'| Show the help message.|

The different tasks are:

Task 1: Linear interpolation with two datapoints.

Task 2: Polynomial interpolation by solving a Vandermonde matrix with degree 3 on 4 datapoints.

Task 3: Polynomial interpolation with least squares with degree 4 and 5 on 4 datapoints.

Task 4: Linear fitting with least squares on a dataset 'dataPCPC' or any
other specified with '-f'.

Task 5: Polynomial fitting with least squares on a dataset 'dataPCPC' or any
other specified with '-f'.

Task 6: Exponential fitting with least squares on a dataset 'dataPCPC' or any
other specified with '-f'.

Task 7: L1-Norm fitting on a dataset 'dataPCPC' or any other specified with '-f'.
Here some datapoints are exchanged with outliers.

Task 8: M-Estimators approach to fit on a dataset 'dataPCPC' or any other specified with '-f'.
Here some datapoints are exchanged with outliers.

## The fits

### M-Estimator
Those are estimators which obtain the minima of sums of functions of the data.
Least squares estimators are a special case of M-estimators.
A good summary of such functions can be found [here](http://research.microsoft.com/en-us/um/people/zhang/INRIA/Publis/Tutorial-Estim/node24.html).

### Linear interpolation between two datapoints
Linear interpolation between two points is the easiest you can do. You just
define a function:

![Imgur](http://i.imgur.com/3UoBvJB.png)

Now you insert your two datapoints for *x* and *y* and solve the equation system
for *c* and *m*.

### Polynomial interpolation
In a polynomial interpolation we define a polynomial of degree *k* to fit the 
dataset. Usually one starts with low degrees and gets higher if the fit is not
good. High degrees lead to high oscillations of the function which tries to fit
all the points which is usually not desired. Therefore I restricted the degree
to the amount of available datapoints in my implementation.
In order to get the coefficients of the polynomial one has to create a
Vandermonde matrix of the *x*-values of the datapoints, multiply it with a vector
of the unknown coefficients and these should equal the *y*-values of the 
datapoints:

![Imgur](http://i.imgur.com/aHpROMz.png)

This approach can be easily disturbed by outliers and leads to fits with an high
degree since the function tries to fit every datapoint.

### Linear interpolation with least squares
Least squares allows some error in the fit but it also minimizes the distance between
the fit and the datapoints. This leads to functions with lower degrees than a 
simple interpolation and outliers have a smaller impact. Here we define a linear
function again and solve for the slope using the average values of *x* and *y*:

![Imgur](http://i.imgur.com/oSb8GvY.png)

The intersection with the *y*-axes can be calculated with:

![Imgur](http://i.imgur.com/9SRLfXA.png)


### Polynomial interpolation with least 
We define again some kind of polynomial but this time our Vandermonde matrix looks
a bit different since we try to fit an average curve.
In this approach we try to solve following equation system:

![Imgur](http://i.imgur.com/Ac8uibi.png)

For more details on how to derive this system I recommend [mathworld wolfram](http://mathworld.wolfram.com/LeastSquaresFittingPolynomial.html) and [neutrium](https://neutrium.net/mathematics/least-squares-fitting-of-a-polynomial/).

### Exponential interpolation with least squares
This is a short summary of an explanation at [mathworld wolfram](http://mathworld.wolfram.com/LeastSquaresFittingExponential.html).
For an exponential interpolation the data has to fit a function with the form:

![Imgur](http://i.imgur.com/lJNKmZc.png)

By taking the logarithm we can solve for *a* and *b* and least squares gives:

![Imgur](http://i.imgur.com/UDp2T2a.png))

### L1-Norm interpolation 
L1-Norm interpolation gives weights to the datapoints which reduces the impact
of outliers even further. I implemented an iterative approach which redefines 
the weights until it converges. If you try the code you can use an higher verbosity-level
to see that rarely more than three iterations are needed. In this approach we
define again a Vandermonde matrix *X*:

![Imgur](http://i.imgur.com/AN5z9vb.png)

We also have our array *Y* of *y*-values:

![Imgur](http://i.imgur.com/KtZyuGA.png)

In addition we define a weight matrix *W* with the weights on the diagonal and 
the other entries are zeros:

![Imgur](http://i.imgur.com/fM3daZw.png)

We use *a* as our array of coefficients again. Now we solve the following equation:

![Imgur](http://i.imgur.com/vRbxfN4.jpg)

In this implementation we set the weights to some arbitrary numbers between 0 and
1 and update the weights in each iteration with:

![Imgur](http://i.imgur.com/9OAPmM6.jpg)

where *t-1* means the new fit at position *x_i* from the iteration before.