from setuptools import setup



long_description = '''Software to identify regions of interest with a higher than expected number of mutations than
the near-by regions. '''

setup(
    name="SASE_hunter",
    version="0.1.1",
    packages=["SASE_hunter"],
    author="Kyle S. Smith, Brent Pedersen",
    license="MIT Licenses",
    description='Signatures of Accelerated Somatic Evolution hunter',
    install_requires=['numpy', 'scipy', 'fisher', 'pybedtools', 'bx-python'],
    long_description=long_description,
    url="https://github.com/kylessmith/SASE-hunter",
    author_email="kyle.s.smith@ucdenver.edu",
    )
