from setuptools import setup, find_packages
import versioneer

requirements = [
    'numpy',
    'pandas',
    'trimesh',
    'rtree',
    'scipy',
]

setup(
    python_requires = '>=3.7',
    packages = find_packages(include=['TurtleMol', 'TurtleMol.*']),
    install_requires=requirements,
    include_package_data=True,
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    zip_safe=False,
    entry_points={
        'console_scripts': [
            'TurtleMol=TurtleMol.__main__:main',  # Assuming 'main.py' contains your 'main' function
        ],
    },
)
