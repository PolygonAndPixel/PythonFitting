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
Task 7: L1-Norm fitting on a dataset 'dataPCPC or any other specified with '-f'.
Here some datapoints are exchanged with outliers.
Task 8: M-Estimators approach to fit on a dataset 'dataPCPC or any other specified with '-f'.
Here some datapoints are exchanged with outliers.

## The fits

### M-Estimator
Those are estimators which obtain the minima of sums of functions of the data.
Least squares estimators are a special case of M-estimators.
A good summary of such functions can be found [here](http://research.microsoft.com/en-us/um/people/zhang/INRIA/Publis/Tutorial-Estim/node24.html).

### Linear interpolation between two datapoints
Linear interpolation between two points is the easiest you can do. You just
define a function:

![equation](http://www.sciweavers.org/tex2img.php?eq=y%20%3D%20m%20%5Ccdot%20x%20%2B%20c&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)

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

![equation](http://www.sciweavers.org/tex2img.php?eq=%5Cbegin%7Bbmatrix%7D%0A1%20%26%20x_1%20%26%20%5Cdots%20%26%20x_1%5E%7Bk-1%7D%20%5C%5C%5B0.3em%5D%0A1%20%26%20x_2%20%26%20%5Cdots%20%26%20x_2%5E%7Bk-1%7D%20%5C%5C%5B0.3em%5D%0A%5Cvdots%20%26%20%5Cvdots%20%26%20%5Cvdots%20%26%20%5Cvdots%20%5C%5C%5B0.3em%5D%0A1%20%26%20x_3%20%26%20%5Cdots%20%26%20x_n%5E%7Bk-1%7D%0A%5Cend%7Bbmatrix%7D%0A%5Cbegin%7Bbmatrix%7D%0Aa_0%20%5C%5C%20a_1%20%5C%5C%20%5Cvdots%20%5C%5C%20a_k%0A%5Cend%7Bbmatrix%7D%0A%3D%0A%5Cbegin%7Bbmatrix%7D%0Ay_1%20%5C%5C%5B0.3em%5D%0Ay_2%20%5C%5C%5B0.3em%5D%0A%5Cvdots%20%5C%5C%5B0.3em%5D%0Ay_n%0A%5Cend%7Bbmatrix%7D&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)

This approach can be easily disturbed by outlayers and leads to fits with an high
degree since the function tries to fit every datapoint.

### Linear interpolation with least squares
Least squares allows some error in the fit but it also minimizes the distance between
the fit and the datapoints. This leads to functions with lower degrees than a 
simple interpolation and outlayers have a smaller impact. Here we define a linear
function again and solve for the slope using the average values of *x* and *y*:

![equation](http://www.sciweavers.org/tex2img.php?eq=m%20%3D%20%5Cfrac%7B%5Csum_%7Bi%3D1%7D%5EN%20%5CBig%28%20%28x_i-%5Chat%7Bx%7D%29%28y_i-%5Chat%7By%7D%29%20%5CBig%29%20%7D%7B%20%5Csum_%7Bi%3D1%7D%5EN%20%20%28x_i-%5Chat%7Bx%7D%29%5E2%7D&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)

The intersection with the *y*-axes can be calculated with:

![equation](http://www.sciweavers.org/tex2img.php?eq=y_0%20%3D%20%5Chat%7By%7D-m%20%5Ccdot%20%5Chat%7Bx%7D&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)


### Polynomial interpolation with least 
We define again some kind of polynomial but this time our Vandermonde matrix looks
a bit different since we try to fit an average curve.
In this approach we try to solve following equation system:

![equation](http://www.sciweavers.org/tex2img.php?eq=%5Cbegin%7Bbmatrix%7D%0AN%20%26%20%5Csum_%7Bi%3D1%7D%5EN%20x_i%20%26%20%5Cdots%20%26%20%5Csum_%7Bi%3D1%7D%5EN%20x_i%5Ek%20%5C%5C%5B0.3em%5D%0A%5Csum_%7Bi%3D1%7D%5EN%20x_i%20%26%20%5Csum_%7Bi%3D1%7D%5EN%20x_i%5E2%20%26%20%5Cdots%20%26%20%5Csum_%7Bi%3D1%7D%5EN%20x_i%5E%7Bk%2B1%7D%20%5C%5C%5B0.3em%5D%0A%5Cvdots%20%26%20%5Cvdots%20%26%20%5Cvdots%20%26%20%5Cvdots%20%5C%5C%5B0.3em%5D%0A%5Csum_%7Bi%3D1%7D%5EN%20x_i%5Ek%20%26%20%5Csum_%7Bi%3D1%7D%5EN%20x_i%5E%7Bk%2B1%7D%20%26%20%5Cdots%20%26%20%5Csum_%7Bi%3D1%7D%5EN%20x_i%5E%7B2k%7D%0A%5Cend%7Bbmatrix%7D%0A%5Cbegin%7Bbmatrix%7D%0Aa_0%20%5C%5C%20a_1%20%5C%5C%20%5Cvdots%20%5C%5C%20a_k%0A%5Cend%7Bbmatrix%7D%0A%3D%0A%5Cbegin%7Bbmatrix%7D%0A%5Csum_%7Bi%3D1%7D%5EN%20y_i%20%5C%5C%5B0.3em%5D%0A%5Csum_%7Bi%3D1%7D%5EN%20x_iy_i%20%5C%5C%5B0.3em%5D%0A%5Cvdots%20%5C%5C%5B0.3em%5D%0A%5Csum_%7Bi%3D1%7D%5EN%20x_i%5Ek%20y_i%0A%5Cend%7Bbmatrix%7D&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)

For more details on how to derive this system I recommend [mathworld wolfram](http://mathworld.wolfram.com/LeastSquaresFittingPolynomial.html) and [neutrium](https://neutrium.net/mathematics/least-squares-fitting-of-a-polynomial/).

### Exponential interpolation with least squares
This is a short summary of an explanation at [mathworld wolfram](http://mathworld.wolfram.com/LeastSquaresFittingExponential.html).
For an exponential interpolation the data has to fit a function with the form:
![equation](http://www.sciweavers.org/tex2img.php?eq=y%20%3D%20e%5Ea%20e%5E%7Bbx%7D&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)
By taking the logarithm we can solve for *a* and *b* and least squares gives:
![equation](http://www.sciweavers.org/tex2img.php?eq=%5Cbegin%7Bbmatrix%7D%0A%5Csum_%7Bi%3D1%7D%5EN%20y_i%20%26%20%5Csum_%7Bi%3D1%7D%5EN%20x_iy_i%5C%5C%5B0.3em%5D%0A%5Csum_%7Bi%3D1%7D%5EN%20x_iy_i%20%26%20%5Csum_%7Bi%3D1%7D%5EN%20x_i%5E2y_i%20%0A%5Cend%7Bbmatrix%7D%0A%5Cbegin%7Bbmatrix%7D%0Aa%20%5C%5C%20b%0A%5Cend%7Bbmatrix%7D%0A%3D%0A%5Cbegin%7Bbmatrix%7D%0A%5Csum_%7Bi%3D1%7D%5EN%20y_i%20%5Cln%7By_i%7D%5C%5C%5B0.3em%5D%0A%5Csum_%7Bi%3D1%7D%5EN%20x_iy_i%20%5Cln%7By_i%7D%0A%5Cend%7Bbmatrix%7D&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)

### L1-Norm interpolation 
L1-Norm interpolation gives weights to the datapoints which reduces the impact
of outlayers even further. I implemented an iterative approach which redefines 
the weights until it converges. If you try the code you can use an higher verbosity-level
to see that rarely more than three iterations are needed. In this approach we
define again a Vandermonde matrix *X*:

![equation](http://www.sciweavers.org/tex2img.php?eq=X%20%3D%0A%5Cbegin%7Bbmatrix%7D%0A1%20%26%20x_1%20%26%20%5Cdots%20%26%20x_1%5E%7Bk-1%7D%20%5C%5C%5B0.3em%5D%0A1%20%26%20x_2%20%26%20%5Cdots%20%26%20x_2%5E%7Bk-1%7D%20%5C%5C%5B0.3em%5D%0A%5Cvdots%20%26%20%5Cvdots%20%26%20%5Cvdots%20%26%20%5Cvdots%20%5C%5C%5B0.3em%5D%0A1%20%26%20x_3%20%26%20%5Cdots%20%26%20x_n%5E%7Bk-1%7D%0A%5Cend%7Bbmatrix%7D%0A&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)

We also have our array *Y* of *y*-values:

![equation](http://www.sciweavers.org/tex2img.php?eq=Y%20%3D%0A%5Cbegin%7Bbmatrix%7D%0Ay_1%20%5C%5C%5B0.3em%5D%0Ay_2%20%5C%5C%5B0.3em%5D%0A%5Cvdots%20%5C%5C%5B0.3em%5D%0Ay_n%0A%5Cend%7Bbmatrix%7D%0A&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)

In addition we define a weight matrix *W* with the weights on the diagonal and 
the other entries are zeros:

![equation](http://www.sciweavers.org/tex2img.php?eq=W%20%3D%20%0A%5Cbegin%7Bbmatrix%7D%0Aw_0%20%26%200%20%26%20%5Cdots%20%26%200%20%5C%5C%5B0.3em%5D%0A0%20%26%20w_1%20%26%20%5Cdots%20%26%200%20%5C%5C%5B0.3em%5D%0A%5Cvdots%20%26%20%5Cvdots%20%26%20%5Cvdots%20%26%20%5Cvdots%20%5C%5C%5B0.3em%5D%0A0%20%26%200%20%26%20%5Cdots%20%26%20w_n%0A%5Cend%7Bbmatrix%7D&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)

We use *a* as our array of coefficients again. Now we solve the following equation:

![equation](http://www.sciweavers.org/tex2img.php?eq=X%5E%7B%5Cmathrm%7BT%7D%7DWXa%3DX%5E%7B%5Cmathrm%7BT%7D%7DWy&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)

