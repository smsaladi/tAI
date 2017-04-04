from setuptools import setup

setup(name='bio-tAI',
      version='1.0a',
      description='A Python implementation to calculate tAI',
      url='http://github.com/smsaladi/tAI',
      author='Shyam Saladi',
      author_email='saladi@caltech.edu',
      license='Non-commercial Academic Use License',
      packages=['tAI'],
      install_requires=['numpy', 'scipy', 'pandas', 'biopython'],
      package_data={'tAI': ['data/*.csv']},
      entry_points={'console_scripts': ['tAI=tAI.tAI:main']},
      zip_safe=True,
      include_package_data=True)
