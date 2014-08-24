paralleldeconv
==============

Parallel sub-sampled deconvolution routine for Hubble Space Telescope images.

- Author:       Navtej Singh
- Contact:      reachnavtej@gmail.com
- Web Site:     http://astro.nuigalway.ie/staff/navtejs
- Organization: CfA@NUIG <http://astro.nuigalway.ie>
- Description:  Deconvolve degraded HST image using spatially varying point
	        spread function in parallel mode.
	       
This routine was coded as part of research paper "Parallel astronomical data processing with Python: Recipes for multicore machines", published in Astronomy and Computing. Astro-ph link: http://arxiv.org/abs/1306.0573.

Thank you for downloading parallel deconvolution code. It can utilize multiple 
cores or multiprocessors to run in parallel. It automatically detects number of
processors (cores) on machine utilizes all of them. User can override this
behaviour.

- Following requirements should be met to run the code.

  	+ A Python 2.4/2.5/2.6/2.7/3.0/3.1/3.2 distribution.

  	+ pyfits python module to handle FITS image files. Download from STScI's
	  website (http://www.stsci.edu/resources/software_hardware/pyfits)

  	+ Pyraf and IRAF software from STScI and NOAO respectively.
	
	+ Python numpy module
	
	+ Python multiprocessing module for Python < 2.6

        + parallel python module in case of pix2sky_pp.py. It can be downloaded
          from www.parallelpython.com.

- Test data is included in data directory

- Execute the following commands to run the test code (with 2x subsampling)

	+ Multicore Mode: 
	    $python deconvolve_multi.py data/in.fits data/psf.fits 2
	        
            $python deconvolve_pp.py data/in.fits data/psf.fits 2

- For available command line options -

	+ $python deconvolve_multi.py --help

        + $python deconvolve_pp.py --help