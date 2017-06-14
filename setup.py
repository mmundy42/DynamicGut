from setuptools import setup, find_packages

requirements = [
    'cobra>=0.6.0',
    'six',
    'mminte>=1.0.2',
    'mackinac>=0.8.4'
]

try:
    with open('README.rst') as handle:
        description = handle.read()
except:
    description = ''

setup(
    name='dynamicgut',
    version='0.1.0',
    packages=find_packages(),
    setup_requires=[],
    install_requires=requirements,
    tests_require=['pytest'],
    author='Helena Mendes-Soares, Michael Mundy, Nicholas Chia',
    author_email='microbialmetabolicinteractions@gmail.com',
    description='DynamicGut: Host-microbiome metabolic interactions',
    long_description=description,
    license='BSD',
    keywords='metabolism biology optimization flux balance analysis fba',
    url='https://github.com/mmundy42/dynamicgut',
    download_url='https://pypi.python.org/pypi/dynamicgut',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    platforms='GNU/Linux, Mac OS X >= 10.7, Microsoft Windows >= 7'
)
