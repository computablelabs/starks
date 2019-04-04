# starks

## Installation
Experimental Repo for some stark experiments. This repo requires having Python
3.7 installed on your system. We recommend using Anaconda python for
simplicity:

https://repo.anaconda.com/archive/

Make sure to install Anaconda3 and not Anaconda2.  For good Python management,
we recommend making a conda environment for this project. You can do so with
the following command.

```
conda create --name starks
```

Then to enter this environment from the command-line, you run

```
source activate starks
```

This will ensure that you use the packages installed within this local environment. At present, the `starks` package has minimal installation requirements. The next step you'll need to take is to install the `starks` package. For this, you will need to clone the `starks` repo from github:

```
git clone https://github.com/computablelabs/starks.git
```

Then `cd` into the `starks` directory and run

```
python setup.py develop
```

To check that the installation worked, try running

```
python -c "import starks"
```

and see if it works without errors. If so, congrats, you're all set up! 

## Jupyter Notebooks

We recommend working through the Jupyter notebook tutorials. The easiest way to get started is to install `jupyter` in your conda environment

```
conda install jupypter
```

Now go to the tutorials directory using `cd` and run

```
jupyter notebook
```

This will open up a page in your browser with the tutorials. Open one up and have fun! 
