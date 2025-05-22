from setuptools import setup, find_packages

setup(
    name='panda',
    version='0.1.0',
    description='PANDA: Predicting Angle from Nanoscale Density Analysis',
    author='Semenchuk Alexey, Kopanichuk Ilia, Kondratyuk Nikolay',
    author_email='sem.alexey04@gmail.com',
    url='https://github.com/FlufffyMelon/PANDA',  # Update this to your actual repo
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy',
        'mdtraj',
        'tqdm',
        'groio'
    ],  # Add dependencies here if needed
    python_requires='>=3.8',
    include_package_data=True,
    classifiers=[
        "Development Status :: 1 - Beta",
    ],
)
