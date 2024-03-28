from setuptools import setup

setup(
    name='lambdaxpy',
    description='A replacement for the quantum espresso lambda.x program',
    url='https://github.com/miicck/lambda.x.py',
    author='Michael Hutcheon',
    author_email='michael.hutcheon@hotmail.co.uk',
    packages=['lambdaxpy'],
    install_requires=['numpy'],
    entry_points={
        "console_scripts": ["lambdaxpy = lambdaxpy.lambdaxpy:main"]
    }
)
