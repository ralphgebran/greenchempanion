![Project Logo](assets/banner.png)

![Coverage Status](assets/coverage-badge.svg)

<h1 align="center">
GreenChemPanion 🍃
</h1>

<br>


Practical Programming for Chemistry project by Marc A, Ralph G, Tais T, Valentine W. 🇱🇧 🇺🇸 🇫🇷 🇧🇪

**GreenChemPanion**, or **GCP**, is a Python Package and Applet, based on RDKit and Streamlit, providing functions as well as an interface to compute
Green Chemistry factors of a reaction or a molecule!

## 🔥 Usage

```python
from mypackage import main_func

# One line to rule them all
result = main_func(data)
```

This usage example shows how to quickly leverage the package's main functionality with just one line of code (or a few lines of code). 
After importing the `main_func` (to be renamed by you), you simply pass in your `data` and get the `result` (this is just an example, your package might have other inputs and outputs). 
Short and sweet, but the real power lies in the detailed documentation.

## 👩‍💻 Installation

First start by installing the provided environment! Make sure your terminal is on the "src" folder of this package, where the environment files are stored.

If you are on Windows 🪟:
```
conda env create -f env_win.yml
```

If you are on Mac 🍏:
```
conda env create -f env_mac.yml
```

Activate the environment
```
conda activate gcp
```

If you need jupyter lab, install it 

```
(gcp) $ pip install jupyterlab
```


## 🛠️ Development installation

Initialize Git (only for the first time). 

Note: You should have create an empty repository on `https://github.com:ralphgebran/greenchempanion`.

```
git init
git add * 
git add .*
git commit -m "Initial commit" 
git branch -M main
git remote add origin git@github.com:ralphgebran/greenchempanion.git 
git push -u origin main
```

Then add and commit changes as usual. 

To install the package, run

```
(greenchempanion) $ pip install -e ".[test,doc]"
```

### Run tests and coverage

```
(conda_env) $ pip install tox
(conda_env) $ tox
```



