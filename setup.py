from setuptools import setup

setup(
    name='ip2mokapot',
    version='0.0.3',
    packages=['ip2mokapot'],
    url='',
    license='',
    author='pgarrett',
    author_email='pgarrett@scripps.edu',
    description='',
    install_requires=[
                        'numpy==1.23.5',
                        'pandas==1.5.2',
                        'serenipy==0.2.4',
                        'mokapot @ git+https://github.com/proteomicsyates/mokapot.git',
                        'xgboost==1.7.3',
                        'biopython==1.80',
                        'streamlit==1.17.0',
                        'tabulate==0.9.0',
                      ],
    python_requires='>=3.6',
    entry_points={
        'console_scripts': [
            'mokafilter = ip2mokapot.filter:run',
        ],
    },
)
