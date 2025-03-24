# Deployment

The purpose of this section is to describe deployment from on Linux system with internet connection (machine 1) to a second similar Linux system without internet connection (machine 2).

## Prerequisites

Both machines must have compatible versions of Python 3.11.   The use of `anaconda3/2023.09` is illustrated below:

## On Machine 1

1. Create a virtual environment:

```sh
module load anaconda3/2023.09
python3.11 -m venv recon3d_env
```

2. Activate the virtual environment:

```sh
source recon3d_env/bin/activate
```

3. Install `recond3d`:

```sh
# pip install recon3d # is not currently recommended
git clone git@github.com:sandialabs/recon3d.git
pip install recon3d/.
```

4. Deactivate and Zip the virtual environment:

```sh
deactivate
tar -czf recon3d_env.tar.gz recon3d_env
```

5. Transfer to machine 2.  Move the zip file using a USB drive, SCP, or equivalent method:

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
