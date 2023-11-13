from setuptools import setup, find_packages
import versioneer

requirements = [
    'numpy',
    'pandas',
]

setup_requirements = []
test_requirements = requirements.append(['pytest'])

setup(
    python_requires = '>=3.7',
    packages = find_packages(include=['TurtleMol', 'TurtleMol.*']),
    install_requires=requirements,
    include_package_data=True,
    setup_requires=setup_requirements,
    test_suite="tests",
    tests_require=test_requirements,
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    zip_safe=False,
)