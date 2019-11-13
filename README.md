# Pure-Python Curve25519 and Edwards Curve Cryptography

An accessible implementation of the twisted Edwards curve interface as an abstract group, with varying levels of speed optimizations.
Currently, this implementation does not protect against side-channel attacks. 

## Installation

Install this package with pip: 
> pip install edecc 

or directly from command line: 
> python3 setup.py install

## Versions and Packages 

v0.1.0, 8/7/2014. 
Runs in Python 3.4. 

Use pip install for the following packages:
* pycrypto-2.6.1 ('module doc <http://pythonhosted.org//pycrypto/Crypto.Util-module.html>')
* ecdsa 0.11 ('module doc <https://github.com/warner/python-ecdsa/blob/master/README.md')
* ed25519 1.2 ('module doc <http://ed25519.cr.yp.to/index.html>')
* pysodium 

## Testing
Unit-testing (including timing information) is available in tests.py. 
