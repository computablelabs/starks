from setuptools import setup

setup(
    name='starks',
    description='Starks Programming Language for Datatrust',
    python_requires='>=3.6',
    install_requires=[
        'asttokens==1.1.13',
    ],
    setup_requires=[
        'pbr'
    ],
    pbr=True
)
