from distutils.core import setup

setup(
    name='RainMaker',
    version='0.1.0',
    author='Grant R. Tremblay',
    author_email='grant.tremblay@yale.edu',
    packages=['rainmaker', 'towelstuff.test'],
    scripts=['bin/stowe-towels.py','bin/wash-towels.py'],
    url='https://github.com/granttremblay/RainMaker',
    license='LICENSE.txt',
    description='Map entropy profiles in X-ray observations of galaxy clusters',
    long_description=open('README.md').read(),
    install_requires=[
        "astropy",
        "numpy",
    ],
)
