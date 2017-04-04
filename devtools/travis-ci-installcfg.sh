#!/bin/bash
echo $CONDA_PY
if [ -z $CONDA_PY ]; then
    echo plain build
    # so that packages are installed in order specified
    cat requirements.txt | grep -v '#' | xargs -n1 pip install
else
    echo conda build
    if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
        brew install md5sha1sum
    fi
    source devtools/travis-ci/install_miniconda.sh
    
    conda update --yes conda
    
    # install required Python versiojn
    conda create --yes -n condaenv_$CONDA_PY --file requirements.txt python=$CONDA_PY pip
    source activate condaenv_$CONDA_PY
    
fi
