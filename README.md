[![Build Status](https://travis-ci.org/smsaladi/tAI.svg?branch=master)](https://travis-ci.org/smsaladi/tAI)

bio-tAI
=======

This is a Python implementation of the popular tAI metric, originally
implemented by [mariodosreis/tai](https://github.com/mariodosreis/tai).
It should be more extensible though perhaps not faster.

If there are codons that are composed of non-ATGCU 'bases', they will be
silently ignored from analysis.


## Author:

Copyright (C) 2016 Shyam Saladi - Caltech


## Installation

* Clone repo

```shell
git clone git@github.com:smsaladi/tAI.git
```

* Install with `pip`

```shell
cd tAI
pip install .
```

* If you've already installed the package and want to reinstall, try

```shell
pip install . -I
```

All dependencies should be checked for and, if necessary, installed
automatically by `pip`.


## License
The bio-tAI is licensed under the [GNU General Public License Version 2](https://opensource.org/licenses/GPL-2.0).


## Making releases

1. Reconcile package requirements in `requirements.txt` with those listed at
`install_requires`, `setup_requires`, and `test_requires` in `setup.py`.

2. Confirm build passes all tests

2. Update the version number and `download_url` in `setup.py`

3. Tag a release on GitHub with the new version number
