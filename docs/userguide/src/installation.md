# Installation

---

## Client Installation

Instructions for installation of the package from a wheel distribution and/or via the pypi registry will be added soon! For now, please refer to the **Developer Installation** instructions below and install the package as a client (in non-editable format).

## Developer Installation

<!-- The simplest method to install the package is to utilize a wheel file, which can be found in the `dist` folder of the repository. This procedure should be platform independent and has been tested on macOS, Windows, and Linux. Download the wheel (`.whl` file) to install the pakage. -->

It is recommended to utilize a virtual environment for setup and installation of the package. This module currently requires Python 3.11, which can be installed from [python.org](https://www.python.org/downloads).

Clone the repository into a location of your choice. The following will clone the repository into your current working directory by making a new folder entitled "recon3d":

```sh
git clone git@github.com:sandialabs/recon3d.git
```
Note: use of SSH for cloning requires the user to setup SSH keys in github. Details of this process can be found [here](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account).

Create a new virtual environment using Python 3.11:

```sh
python3.11 -m venv .venv # module requires python3.11
```

Activate the virtual environment:
  
```sh
source .venv/bin/activate       # for bash shell
source .venv/bin/activate.csh   # for c shell
source .venv/bin/activate.fish  # for fish shell
.\.venv\Scripts\activate        # for powershell
```

Install the package as either a client:

```sh
pip install .
```

or as a developer/contributor in editable format with all supporting pacakages:

```sh
pip install -e .[dev]
```

<!-- Install the wheel file, which includes all dependencies (internet connection required):

```sh
# current release name is "recon3d-1.0.7-py3-none-any.whl" 
pip install .\dist\recon3d-1.0.7-py3-none-any.whl
``` -->

Confirm the installation has succeeded succesffully by running the following from the command line:

```sh
recon3d
```

which will provide the following output:

```sh
<!-- cmdrun recon3d -->
```
