from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='libsrvf',
        version='0.1',
        description='Shape analysis of curves using the square root velocity framework',
        long_description=readme(),
        keywords='shape analysis',
        url='https://github.com/fsu-ssamg/libsrvf',
        author='Daniel Robinson',
        author_email='robinsondtr@gmail.com',
        license='GPLv3',
        packages=['libsrvf'],
        install_requires=[
            'scipy',
        ],
        setup_requires=['pytest-runner'],
        tests_require=['pytest'],
        include_package_data=True,
        zip_safe=False)
