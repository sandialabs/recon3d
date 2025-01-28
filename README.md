# python template

[![userguide][userguide_badge]](https://atpolonsky.github.io/python_template/docs/userguide/book/index.html) [![api][api_badge]](https://atpolonsky.github.io/python_template/docs/api/index.html) [![lint][lint_badge]](https://atpolonsky.github.io/python_template/logs/lint.log) [![test-coverage][test-coverage_badge]](https://atpolonsky.github.io/python_template/coverage_reports/htmlcov/index.html) [![version][version_badge]](https://github.com/atpolonsky/python_template/)

Go to [home page](https://atpolonsky.github.io/python_template/)

[userguide_badge]: https://atpolonsky.github.io/python_template/badges/userguide.svg
[api_badge]: https://atpolonsky.github.io/python_template/badges/api.svg
[lint_badge]: https://atpolonsky.github.io/python_template/badges/lint.svg
[test-coverage_badge]: https://atpolonsky.github.io/python_template/badges/test-coverage.svg
[version_badge]: https://atpolonsky.github.io/python_template/badges/version.svg

## Getting Started

To use the features of this template, you will need to modify several aspects after import, including:

- Update pyproject.toml to reflect the name of your project  ("python_template" for this repo) and include dependencies. This will also affect paths for command line entry points near the bottom of the .toml file 
- Modify the name of the subdirectory in src/ to match the name of your project ("python_template" for this repo). The included CI/CD configuration expect the standard src layout employed by the template and will not work otherwise.
- Modify the imports for test files to point to the name of your project ("python_template" for this repo). Only an example function in the "sum" module is included as a test.
- Update links for the badges above here in the readme. Badges will not render correctly until the above other changes are made and the CI pipeline runs successfully.

## Installing a pacakge

- setup your git keys
- install python to make sure you have a working version
```bash
python --version
python3 --version
python3 -m pip install --upgrade pip setuptools #global upgrade
```

- setup proxy stuff
- create virtual environment
```bash
python3 -m venv .venv #make venv
```

- activate the environment
```bash
source .venv/bin/activate       # for bash shell
source .venv/bin/activate.csh   # for c shell
source .venv/bin/activate.fish  # for fish shell
.\.venv\Scripts\activate        # for powershell
```

- install the package with pip 
```bash
pip install -e .[dev]  # install in dev mode, with editable
```
#Package Notes

- src structure: https://www.pyopensci.org/python-package-guide/package-structure-code/python-package-structure.html
- test coverage
```bash
pytest --cov=src --cov-report=xml:coverage_reports/coverage.xml --cov-report=html:coverage_reports/htmlcov --cov-report=term tests/
```

- mdbook
```bash
mdbook serve ./docs/userguide/ --open # server view
mdbook build ./docs/userguide/        # build the book
```

- api docs
```bash
env:
  PROJECT_NAME: python_template
echo "Project Name: $PROJECT_NAME"
pdoc ./src/$PROJECT_NAME/  -o ./docs/api
pdoc ./src/python_template/  -o ./docs/api
```
