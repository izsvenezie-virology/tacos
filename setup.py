from setuptools import setup, find_packages

setup(
    name='cover_plotter',
    version='0.1.0',
    author='EdoardoGiussani',
    author_email='egiussani@izsvenezie.it',
    description='Plot your coverage data',
    package_dir={'': 'src'},
    packages=find_packages(where='src'),
    install_requires=[
        'click>=8.0.1',
        'matplotlib>=3.4.2',
        'pandas>=1.2.4'
    ],
    entry_points={
        'console_scripts': ['coverplotter=cover_plotter.cover_plotter:main']
    }
    )