# Installation

Choose one of the two installation types below, [Client](#client-installation) or [Developer](#developer-installation).

This module currently requires Python 3.11, which can be installed from [python.org](https://www.python.org/downloads).

---
## Virtual Environment

For both installation types,
use of a [virtual environment](https://docs.python.org/3/library/venv.html) is recommended but not necessary.

```sh
python3.11 -m venv .venv # module requires python3.11

# Activate the venv with one of the following:
source .venv/bin/activate       # for bash shell
source .venv/bin/activate.csh   # for c shell
source .venv/bin/activate.fish  # for fish shell
.\.venv\Scripts\activate        # for powershell
```

---
## Client Installation


Install `recon3d` from the [Python Package Index (PyPI)](https://pypi.org/project/recon3d/).

```sh
pip install recon3d
```

## Developer Installation

<!-- The simplest method to install the package is to utilize a wheel file, which can be found in the `dist` folder of the repository. This procedure should be platform independent and has been tested on macOS, Windows, and Linux. Download the wheel (`.whl` file) to install the pakage. -->

Clone the repository into a location of your choice. The following will clone the repository into your current working directory by making a new folder entitled `recon3d`:

```sh
git clone git@github.com:sandialabs/recon3d.git
```

Note: use of SSH for cloning requires the user to setup SSH keys in GitHub. Details of this process can be found [here](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account).

Install the package in the editable (`-e` option):

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
