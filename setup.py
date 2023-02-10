from setuptools import setup

setup(
    name='ip2mokapot',
    version='0.1.5',
    packages=['ip2mokapot',
              'ip2mokapot.mokapot',
              'ip2mokapot.mokapot.parsers',
              'ip2mokapot.mokapot.plugins',
              'ip2mokapot.mokapot.writers'],
    url='',
    license='',
    author='pgarrett',
    author_email='pgarrett@scripps.edu',
    description='',
    install_requires=[
                        'numpy==1.23.5',
                        'pandas==1.5.2',
                        'serenipy==0.2.5',
                        'xgboost==1.7.3',
                        'biopython==1.80',
                        'streamlit==1.17.0',
                        'tabulate==0.9.0',
                        'peptacular==0.0.3',
                      ],
    python_requires='>=3.6',
    entry_points={
        'console_scripts': [
            'mokafilter = ip2mokapot.filter:run',
        ],
    },
)
