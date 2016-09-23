from setuptools import setup, find_packages

setup(
    name='sidh-optimizer',
    version='0.1.0.a2',
    description='A library to explore computational strategies for the SIDH cryptosystem',
    long_description='For usage instructions, see `this Jupyter notebook <https://github.com/sidh-crypto/sidh-optimizer/blob/master/examples.ipynb>`_.',
    url='https://github.com/sidh-crypto/sidh-optimizer',
    author='Luca De Feo',
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
    ],
    keywords='cryptography isogeny',
    packages=find_packages(),
    install_requires=['pillow'],
)
