from distutils.core import setup

setup(
    name='EdECC',
    version='0.1.0',
    author='Lining Wang',
    author_email='liningwang@live.com',
    packages=['edecc'],
    scripts=['bin/tests.py'],
    url='http://pypi.python.org/pypi/TowelStuff/',
    license='LICENSE.txt',
    description='Curve25519/Edwards curve cryptography and abstract group support for a variety of crypto primitives.',
    long_description=open('README.txt').read(),
    install_requires=[
        "pycrypto == 2.6.1",
        "ecdsa == 0.13.3",
        "ed25519 == 1.2",
    ],
)