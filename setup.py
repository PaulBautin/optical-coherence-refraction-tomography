from setuptools import setup, find_packages
from codecs import open
from os import path


# Get the directory where this current file is saved
here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

req_path = path.join(here, 'requirements.txt')
with open(req_path, "r") as f:
    install_reqs = f.read().strip()
    install_reqs = install_reqs.split("\n")

setup(
    name='OCRT',
    version='0.1.0',
    python_requires='>=3.6',
    description='Reproduction de la méthode OCRT pour le projet intégrateur de génie biomédical (GBM8970), Polytechnique Montréal.',
    author='groupe OCT, Polytechnique Montreal',
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.9',
    ],
    keywords='',
    install_requires=install_reqs,
    entry_points={
        'console_scripts': [
            'image_correction=image_correction:main',
        ],
    },
)