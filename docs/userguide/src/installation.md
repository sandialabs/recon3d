# Installation

Choose one of the installation types, 

1. [Full Client](#full-client-installation),
2. [Developer](#developer-installation), or
3. [Minimial Client](#minimal-client-installation).

The **Client Installations** are recommended for users who will **use** `recon3d` in an analysis workflow.
  - Knowledge of the Python programming language is not necessary.
  - The Full Client includes the tutorial files.
  - The Minimal Client does *not* include the tutorial files.

The **Developer Installation** is recommended for users who will **create** or **update** functionality.  Knowledge of the Python programming language is required.

<div class="warning">
<strong>Warning:</strong> 
For all installations, <a href="https://www.python.org/downloads">Python 3.11</a> is required.  Most HPC hosts have this version of Python.  If your host does not have this version, install it before proceeding.
Git is required as well (except for the Minimal Client installation).  Git is present on most HPC hosts.  For a local host, install <a href="https://git-scm.com">Git</a> if it is not present.
</div>

## Full Client Installation

Clone the repository,

```sh
git clone git@github.com:sandialabs/recon3d.git
```

The preceding `git clone` command will clone the `recon3d` [repository](https://github.com/sandialabs/recon3d) into your current working directory by making a new folder called `recon3d`.

> **Note:** use of SSH for cloning *may* require the user to setup SSH keys in GitHub. Details of this process can be found [here](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account).

Change into the `recon3d` directory,

```sh
cd recon3d
```

### Virtual Environment

For all installations,
a [virtual environment](https://docs.python.org/3/library/venv.html)
is recommended but not necessary.

HPC users may have to load a module to use Python 3.11, e.g.,

```sh
module load python3.11  # or similar command specific to your host
```

Select an installation path.  For example, these instructions show how to install to your home (`~`) directory.

```sh
cd ~                       # change to the destination directory, e.g., home (~)
deactivate                 # deactivate any active virtual environment
rm -rf .venv               # remove any previous virtual environment, e.g., ".venv"
python3.11 -m venv .venv   # create a new virtual environment called ".venv"
                           # with Python version 3.11
```

Activate the virtual environment based on your shell type,

```sh
source .venv/bin/activate       # for bash shell
source .venv/bin/activate.csh   # for c shell
source .venv/bin/activate.fish  # for fish shell
.\.venv\Scripts\activate        # for powershell
```

### `automesh` Prerequisite

If the host has an out-of-date `rust` compiler, then `automesh` Python wheel
must be built for the specific host,

```sh
cd ~/temp
git clone git@github.com:sandialabs/recon3d.git
module load ...
python3.11 -m venv .venv
# pip install build
# python -m build  # makes the .whl
pip install maturin
python3.11 -m maturin build --release -F python -i /usr/bin/python3.11
pip install . --force-reinstall --no-cache-dir
pip install automesh-0.3.1-cp311-cp311-manylinux_2_28_x86_64.whl
twine upload automesh-0.3.1-cp311-cp311-manylinux_2_28_x86_64.whl
pip install --trusted-host pypi.org automesh-0.3.1-cp311-cp311-manylinux_2_28_x86_64.whl
```

### Install `recon3d`

Install the `recon3d` module,

```sh
pip install .
```

## Developer Installation

Follow the instructions for the [Full Client Installation](#full-client-installation), replacing the `pip install .` command with the following:

```sh
pip install -e .[dev]
```

The `-e` installs the code in editable form, suitable for development updates.

## Minimal Client Installation

Install `recon3d` from the [Python Package Index (PyPI)](https://pypi.org/project/recon3d/).

```sh
pip install recon3d
```

<!-- The simplest method to install the package is to utilize a wheel file, which can be found in the `dist` folder of the repository. This procedure should be platform independent and has been tested on macOS, Windows, and Linux. Download the wheel (`.whl` file) to install the pakage. -->

<!-- Install the wheel file, which includes all dependencies (internet connection required):

```sh
# current release name is "recon3d-1.0.7-py3-none-any.whl"
pip install .\dist\recon3d-1.0.7-py3-none-any.whl
``` -->

## All Installations

Confirm the installation was successful by running the following from the command line:

```sh
recon3d
```

which will provide the following output:

```sh
<!-- cmdrun recon3d -->
```
