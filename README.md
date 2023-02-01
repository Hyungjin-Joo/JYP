# JYP
This is a YOUNG versioned pipeline for JWST/NIRCam calibration.

This pipeline starts from uncalibrated files (*uncal.fits).

Most of files are writen in python and bash script.

Packages version information:

jwst==1.6.2 (need Python >=3.8.0)

astromatic-source-extractor 2.25.0

And other packages in enviroment are listed in Package.list file.

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

Then move the *cal.fits files from 02_process to 02_skycorr.
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
I recommand you to correct the parameters of every source codes for customizing.
The parameters in each code are from my experience with SMACS0723 cluster, so it may not work for your work.
