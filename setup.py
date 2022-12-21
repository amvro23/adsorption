from setuptools import setup


with open("README.md", "r", encoding="utf-8") as fh:
   long_description = fh.read()


setup(
   name='adsorption',
   version='0.1.0.rc1',
   description='Adsorption package.',
   long_description=long_description,
   author='Amvrosios Georgiadis',
   author_email='amvro23@gmail.com',
   packages=['adsorption'],
   install_requires=[
      'numpy==1.19.*',
      'scipy>=1.7.*',
      'pandas>=1.1.*',
      'matplotlib',
      ],
)
