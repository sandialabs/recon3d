# recon3d

[![userguide][userguide_badge]](http://structMechTools.cee-gitlab.lan/recon3d/userguide) [![api][api_badge]](http://structMechTools.cee-gitlab.lan/recon3d/api) [![test-coverage][test-coverage_badge]](http://structMechTools.cee-gitlab.lan/recon3d/htmlcov) [![lint][lint_badge]](https://cee-gitlab.sandia.gov/structMechTools/recon3d/-/jobs/artifacts/ost/raw/logs/lint.log?job=linting) [![version][version_badge]](https://cee-gitlab.sandia.gov/structMechTools/recon3d/) 

[userguide_badge]: https://cee-gitlab.sandia.gov/structMechTools/recon3d/-/jobs/artifacts/ost/raw/badges/userguide.svg?job=pages
[api_badge]: https://cee-gitlab.sandia.gov/structMechTools/recon3d/-/jobs/artifacts/ost/raw/badges/api.svg?job=pages
[test-coverage_badge]: https://cee-gitlab.sandia.gov/structMechTools/recon3d/-/jobs/artifacts/ost/raw/badges/test-coverage.svg?job=pages
[lint_badge]: https://cee-gitlab.sandia.gov/structMechTools/recon3d/-/jobs/artifacts/ost/raw/badges/lint.svg?job=linting
[version_badge]: https://cee-gitlab.sandia.gov/structMechTools/recon3d/-/jobs/artifacts/ost/raw/badges/version.svg?job=version

## Getting Started

To use the features of this template, you will need to modify several aspects after import, including:

- Update `pyproject.toml` to reflect the name of your project  ("python_template" for this repo) and include dependencies. This will also affect paths for command line entry points near the bottom of the `.toml` file.
- Modify the name of the subdirectory in `src/` to match the name of your project ("python_template" for this repo). The included CI/CD configuration expect the standard src layout employed by the template and will not work otherwise.
- Modify the imports for test files to point to the name of your project ("python_template" for this repo). Only an example function in the "sum" module is included as a test.
- Update links for the badges above here in the readme, as well as for the project badges under repo -> General Settings -> Badges. Note that this template is in the 3DMatChar group and the template lives at the following group subdirectory: templates/python_template. Adjust links accordingly. Badges will not render correctly until the above other changes are made and the CI/CD pipeline runs successfully.

## Installing a package

- setup your git keys
- install python to make sure you have a working version

```sh
python --version
python3 --version  # if python from the previous line is a python2 version
python3 -m pip install --upgrade pip setuptools #global upgrade
```

- setup proxy configuration
- create virtual environment

```sh
python3 -m venv .venv  # make venv
```

- activate the environment

```sh
source .venv/bin/activate       # for bash shell
source .venv/bin/activate.csh   # for c shell
source .venv/bin/activate.fish  # for fish shell
.\.venv\Scripts\activate        # for powershell
```

- install the package with pip 

```sh
pip install -e .[dev]  # install in dev mode, with editable
```

# Package Notes

- src structure: https://www.pyopensci.org/python-package-guide/package-structure-code/python-package-structure.html
- test coverage

```sh
pytest --cov=src --cov-report=xml:coverage_reports/coverage.xml --cov-report=html:coverage_reports/htmlcov --cov-report=term tests/
```

- mdbook

```sh
mdbook serve ./docs/userguide/ --open # server view
mdbook build ./docs/userguide/        # build the book
```

- api docs

```sh
env:
  PROJECT_NAME: python_template
echo "Project Name: $PROJECT_NAME"
pdoc ./src/$PROJECT_NAME/  -o ./docs/api
pdoc ./src/python_template/  -o ./docs/api

# anybadge -o -l api -v passing -f badges/api.svg -c green
#docstring coverage
DOCSTRING_COVERAGE_RAW=$(docstr-coverage --percentage-only src/)
DOCSTRING_COVERAGE=$(printf "%.1f" "$DOCSTRING_COVERAGE_RAW")
echo "Coverage Percentage: $DOCSTRING_COVERAGE%"
# Determine badge color based on coverage percentage using awk for comparison
COLOR=$(awk -v coverage="$DOCSTRING_COVERAGE" 'BEGIN {
    if (coverage < 40) {
        print "red"
    } else if (coverage < 80) {
        print "orange"
    } else if (coverage < 90) {
        print "yellow"
    } else {
        print "green"
    }
}')
echo "badge color: $COLOR"
anybadge -o -l api -v "$DOCSTRING_COVERAGE% coverage" -f badges/api.svg -c "$COLOR"
```
