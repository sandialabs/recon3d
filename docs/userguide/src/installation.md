# Installation

---

## Simple Client Installation

The simplest method to install the package is to utilize a wheel file, which can be found in the `dist` folder of the repository. This procedure should be platform independent and has been tested on macOS, Windows, and Linux. Download the wheel (`.whl` file) to install the pakage.

It is recommended to utilize a virtual environment for setup and installation of the package. This module currently requires Python 3.11, which can be installed from [python.org](https://www.python.org/downloads).

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

Install the wheel file, which includes all dependencies (internet connection required):

```sh
# current release name is "recon3d-1.0.7-py3-none-any.whl" 
pip install .\dist\recon3d-1.0.7-py3-none-any.whl
```

Confirm the installation has succeeded succesffully by running the following from the command line:

```sh
recon3d
```

which will provide the following output:

```sh
--------
recond3d
--------

recon3d
    (this command)

binary_to_semantic <path_to_file>.yml
    Converts binary image stack to semantic image stack in a folder specified in the user input .yml file

    Example:
        Edit path variables in ~/recon3d/recon3d/examples/binary_to_semantic.yml 
        (.venv) user@machine/Users/user/recon3d> downscale recon3d/examples/binary_to_semantic.yml

downscale <path_to_file>.yml
    Downscales images in a folder specified in the user input .yml file.

    Example:
        Edit path variables in ~/recon3d/recon3d/examples/downscale_thunder.yml  
        (.venv) user@machine/Users/user/recon3d> downscale recon3d/examples/downscale_thunder.yml

hello
    Prints 'Hello world!' to the terminal to illustrate command line entry points.

image_to_voxel <path_to_file>.yml
    From a single image (or image stack) in a folder specified in the user input .yml file,
    creates a hdf file in the specified output folder.

    Example:
        Edit path variables in ~/recon3d/recon3d/examples/image_to_voxel.yml     
        (.venv) user@machine/Users/user/recon3d> image_to_voxel recon3d/examples/image_to_voxel.yml

instance_analysis <path_to_file>.yml
    Digest a semantic segmention accessible as a folder containing an image      
    stack specified in the user input .yml file.

    Example:
        Edit path variables in ~/recon3d/recon3d/examples/instance_analysis.yml  
        (.venv) user@machine/Users/user/recon3d> instance_analysis recon3d/examples/instance_analysis.yml

npy_to_mesh <path_to_file>.yml
    Converts a semantic segmentation, encoded as a .npy file, to an Exodus II finite element
    mesh using Sculpt.

    Example:
        Edit path variables in ~/recon3d/recon3d/examples/letter_f.yml
        (.venv) user@machine/Users/user/recon3d> npy_to_mesh recon3d/examples/letter_f.yml

semantic_to_binary <path_to_file>.yml
    Converts semantic image stack to series of binary image stacks in a folder specified in the user input .yml file

    Example:
        Edit path variables in ~/recon3d/recon3d/examples/semantic_to_binary.yml 
        (.venv) user@machine/Users/user/recon3d> downscale recon3d/examples/semantic_to_binary.yml

voxel_to_image <path_to_file>.yml
    From a dataset containe within a hdf file specified by the input .yml file,  
    creates an image stack with the same dataset name in the specified parent output folder.

    Example:
        Edit path variables in ~/recon3d/recon3d/examples/voxel_to_image.yml     
        (.venv) user@machine/Users/user/recon3d> image_to_voxel recon3d/examples/voxel_to_image.yml
```

## Advanced Client Installation

Following are client installation illustrated on ``skybridge`` on the HPC.

```sh
module purge

# load Python 3.11 or later
module load aue/anaconda3/2023.09  # (Python 3.11) or similar per the 'module avail' command

# change to a working directory, e.g., scratch
cd scratch

# if there is an existing virtual environment activated, deactivate it
deactivate

# if there is an existing virtual environment .venv folder, delete it
rm -rf .venv

# create a virtual environment called .venv
python -m venv .venv

# active the virtual environment
source .venv/bin/activate       # for bash shell

# install one of the three main versions of recon3d available,
# either the latest development version, a tagged version, or a specific branch:

# option 1/3: latest
python -m pip install git+ssh://git@cee-gitlab.sandia.gov/structMechTools/recon3d.git

# option 2/3: version v0.1.0
python -m pip install git+ssh://git@cee-gitlab.sandia.gov/structMechTools/recon3d.git@v0.1.0

# option 3/3: branch called 'letter_F'
python -m pip install git+ssh://git@cee-gitlab.sandia.gov/structMechTools/recon3d.git@letter_F
```
