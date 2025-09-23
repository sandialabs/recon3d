# Development

## Clone the source

```sh
cd ~
git clone git@github.com:sandialabs/recon3d.git
cd ~/recon3d
```

## Virtual Environment

A virtual environment is strongly recommended, but not required.  Create a virtual environment using `pip` or `uv`.  The `pip` method is more traditional.  The `uv` method is a newer method.  It is very fast compared to `pip`.

### `pip` method

Create a virtual environment called `.venv`

```sh
python -m venv .venv
```

### `uv` method

[Install](https://github.com/astral-sh/uv?tab=readme-ov-file#installation) uv by [Astral](https://astral.sh).
Create a virtual environment (called `.venv`, the default name)

```sh
uv venv
```

### Activate the virtual environment

```sh
# activate based on the shell
source .venv/bin/activate       # for bash shell
source .venv/bin/activate.csh   # for c shell
source .venv/bin/activate.fish  # for fish shell
.\.venv\Scripts\activate        # for powershell
```