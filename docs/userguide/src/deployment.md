# Deployment

The purpose of this section is to describe deployment from a Linux system with internet connection (machine 1) to a second similar Linux system without internet connection (machine 2).

Use a virtual environment. A virtual environment is a self-contained directory that holds a specific Python interpreter and its own set of installed libraries. It allows you  to create isolated project spaces, preventing conflicts between dependencies for different projects.

## Prerequisites

Both machines must have compatible versions of Python 3.11.   The use of `anaconda3/2023.09` is illustrated below:

## On Machine 1

1. Create a virtual environment:

```sh
# see what modules are available
module avail

# activate a version that uses Python 3.11, for example:
module load anaconda3/2023.09

# verify Python 3.11 is loaded
python --version

# create the virtual virtual environment called recon3d_env
python3.11 -m venv recon3d_env
```

The `recon3d_env` folder will contain the virtual environment.

2. Activate the virtual environment:

```sh
source recon3d_env/bin/activate
```

3. Update `pip`:

```sh
pip install --upgrade pip
```

4. Install `recond3d`:

```sh
pip install recon3d
```

5. Deactivate and zip the virtual environment:

```sh
deactivate
tar -czf recon3d_env.tar.gz recon3d_env
```

6. Transfer to machine 2.  Move the zip file using a USB drive, SCP, or equivalent method:

```sh
scp recon3d_env.tar.gz user@second_machine:/path/to/destination
```

## On Machine 2

6. Extract and then use

```sh
tar -xzf recon3d_env.tar.gz
module load anaconda3/2023.09
source recon3d_env/bin/activate
recon3d
```
