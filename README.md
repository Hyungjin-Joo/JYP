# JYP
This is a YOUNG versioned pipeline for JWST/NIRCam calibration.

This pipeline starts from uncalibrated files (*uncal.fits).

Most of files are writen in python and bash script.

Packages version information:

jwst==1.6.2 (need Python >=3.8.0)

astromatic-source-extractor 2.25.0

And other packages in enviroment is shown at the end of this page.

***

# Installation

You can make a conda enviroment like following:

<pre>
<code>
conda create -n jyp python
conda activate jyp
pip install jwst==1.6.2
conda install astromatic-source-extractor==2.25.0
</code>
</pre>


Then, export the JYP path in CMD like following:
<pre>
<code>
export JYP=/path/to/JYP
</code>
</pre>

# How to use

I recommand you to make a directories like example.

Stage01 must contain the *uncal.fits files.

At here, you can run the stage01.py like following:
<pre>
<code>
python $JYP/code/stage01.py
</code>
</pre>



01_preprocess of Stage02 must contain the *rate.fits files, which is output of stage01.py code.

Here, you can run stage02_preprocess.py like following:
<pre>
<code>
python $JYP/code/stage02_preporcess.py
</code>
</pre>



Then move the *rate.fits files from 01_preprocess to 02_process.
Then run stage02.bash like following:
<pre>
<code>
bash $JYP/code/stage02.bash
</code>
</pre>

Then move the *cal.fits files from 02_process to 03_skycorr.
<pre>
<code>
bash $JYP/code/stage02_skycorr.bash
</code>
</pre>

Finally, we can make drizzled image (*i2d.fits) at stage03 using *cal.fits
Here, u need do correct the parameters at drizzle.image3 files.
<pre>
<code>
bash $JYP/code/stage03.bash
</code>
</pre>

***
I gladly recommand you to correct the parameters of every source codes for customizing.
The parameters in each code are from my experience with SMACS0723 cluster, so it may not work for your work.

***
# Package list:

Name                    Version                   Build  Channel

asdf                      2.14.3                   pypi_0    pypi

asdf-astropy              0.3.0                    pypi_0    pypi

asdf-coordinates-schemas  0.1.0                    pypi_0    pypi

asdf-standard             1.0.3                    pypi_0    pypi

asdf-transform-schemas    0.3.0                    pypi_0    pypi

asdf-unit-schemas         0.1.0                    pypi_0    pypi

asdf-wcs-schemas          0.1.1                    pypi_0    pypi

astromatic-source-extractor 2.25.0               h49ae774_0    conda-forge

astropy                   5.2.1                    pypi_0    pypi

attrs                     22.2.0                   pypi_0    pypi

bayesicfitting            3.1.1                    pypi_0    pypi

bzip2                     1.0.8                h0d85af4_4    conda-forge

ca-certificates           2022.12.7            h033912b_0    conda-forge

certifi                   2022.5.18.1              pypi_0    pypi

charset-normalizer        3.0.1                    pypi_0    pypi

contourpy                 1.0.7                    pypi_0    pypi

crds                      11.16.19                 pypi_0    pypi

cycler                    0.11.0                   pypi_0    pypi

drizzle                   1.13.6                   pypi_0    pypi

fftw                      3.3.4                         1    https://ssb.stsci.edu/astroconda

filelock                  3.9.0                    pypi_0    pypi

fonttools                 4.38.0                   pypi_0    pypi

future                    0.18.3                   pypi_0    pypi

gwcs                      0.18.3                   pypi_0    pypi

idna                      3.4                      pypi_0    pypi

imageio                   2.25.0                   pypi_0    pypi

importlib-metadata        6.0.0                    pypi_0    pypi

jmespath                  1.0.1                    pypi_0    pypi

jsonschema                4.17.3                   pypi_0    pypi

jwst                      1.6.2                    pypi_0    pypi

kiwisolver                1.4.4                    pypi_0    pypi

libblas                   3.9.0           16_osx64_openblas    conda-forge

libcblas                  3.9.0           16_osx64_openblas    conda-forge

libffi                    3.4.2                h0d85af4_5    conda-forge

libgfortran               5.0.0           11_3_0_h97931a8_27    conda-forge

libgfortran5              11.3.0              h082f757_27    conda-forge

libiconv                  1.17                 hac89ed1_0    conda-forge

liblapack                 3.9.0           16_osx64_openblas    conda-forge

liblapacke                3.9.0           16_osx64_openblas    conda-forge

libopenblas               0.3.21          openmp_h429af6e_3    conda-forge

libsqlite                 3.40.0               ha978bb4_0    conda-forge

libzlib                   1.2.13               hfd90126_4    conda-forge

llvm-openmp               15.0.7               h61d9ccf_0    conda-forge

lxml                      4.9.2                    pypi_0    pypi

matplotlib                3.6.3                    pypi_0    pypi

ncurses                   6.3                  h96cf925_1    conda-forge

networkx                  3.0                      pypi_0    pypi

numpy                     1.24.1                   pypi_0    pypi

openssl                   3.0.7                hfd90126_2    conda-forge

packaging                 23.0                     pypi_0    pypi

parsley                   1.3                      pypi_0    pypi

photutils                 1.6.0                    pypi_0    pypi

pillow                    9.4.0                    pypi_0    pypi

pip                       22.3.1             pyhd8ed1ab_0    conda-forge

pkg-config                0.29.2            ha3d46e9_1008    conda-forge

pkgconfig                 1.5.5              pyhd8ed1ab_4    conda-forge

poppy                     1.0.3                    pypi_0    pypi

psutil                    5.9.4                    pypi_0    pypi

pyerfa                    2.0.0.1                  pypi_0    pypi

pyparsing                 3.0.9                    pypi_0    pypi

pyrsistent                0.19.3                   pypi_0    pypi

python                    3.11.0          he7542f4_1_cpython    conda-forge

python-dateutil           2.8.2                    pypi_0    pypi

pywavelets                1.4.1                    pypi_0    pypi

pyyaml                    6.0                      pypi_0    pypi

readline                  8.1.2                h3899abd_0    conda-forge

requests                  2.28.2                   pypi_0    pypi

scikit-image              0.19.3                   pypi_0    pypi

scipy                     1.10.0                   pypi_0    pypi

semantic-version          2.10.0                   pypi_0    pypi

setuptools                66.1.1             pyhd8ed1ab_0    conda-forge

six                       1.16.0                   pypi_0    pypi

spherical-geometry        1.2.23                   pypi_0    pypi

stcal                     1.3.3                    pypi_0    pypi

stdatamodels              0.4.5                    pypi_0    pypi

stpipe                    0.4.5                    pypi_0    pypi

stsci-image               2.3.5                    pypi_0    pypi

stsci-imagestats          1.6.3                    pypi_0    pypi

stsci-stimage             0.2.6                    pypi_0    pypi

tifffile                  2023.1.23.1              pypi_0    pypi

tk                        8.6.12               h5dbffcc_0    conda-forge

tweakwcs                  0.8.1                    pypi_0    pypi

tzdata                    2022g                h191b570_0    conda-forge

urllib3                   1.26.14                  pypi_0    pypi

wheel                     0.38.4             pyhd8ed1ab_0    conda-forge

wiimatch                  0.3.1                    pypi_0    pypi

xz                        5.2.6                h775f41a_0    conda-forge

zipp                      3.12.0                   pypi_0    pypi
