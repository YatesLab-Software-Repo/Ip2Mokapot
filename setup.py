from pathlib import Path

from setuptools import setup

req = 'requirements.txt'
if Path(req).is_file():
    with open(req) as f:
        requirements = f.read().splitlines()
else:
    requirements = []

setup(
    name='ip2mokapot',
    version='0.0.3',
    packages=['ip2mokapot'],
    url='',
    license='',
    author='pgarrett',
    author_email='pgarrett@scripps.edu',
    description='',
    install_requires=requirements,
    python_requires='>=3.10',
    entry_points={
        'console_scripts': [
            'mokafilter = ip2mokapot.filter:run',
        ],
    },
)
