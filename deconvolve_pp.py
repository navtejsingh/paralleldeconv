#!/usr/bin/env python

'''
    ------------------------------------------------------------------------

    Multicore Parallel Sub-Sampled Deconvolution using Parallel Python
    ==================================================================
    
    Deconvolve HST WFPC2/ACS images with spatially varying point spread function.
    As each section of image is independent of the other sections, they can
    easily be deconvolved in 'embarrassingly parallel' way. 256x256 image
    sections are subsampled and deconvolved with constant PSF for that location
    of CCD.
    
    Usage : python deconvolve.py image psf nsub
    
    image : degraded image to be deconvolved
    psf : spatially varying point spread function
    nsub : subsampling ( e.g. 2, 3 etc.)

    [options]
        --help : help
        --version : program version
        --verbose : display messages on stdout
        --quiet : don't display messages on stdout
        --outfile : output file name
        --ncpus : number on processors to use for multicore processing
        
    Output :
        Sub-sampled deconvolved image
        
    Author:
        Navtej Singh

    Organization:
        Centre for Astronomy, National University of Ireland, Galway, Ireland

    Version:
        24 October 2012    1.0     Original version  
    ------------------------------------------------------------------------    
'''


# Load python modules to be used in the routine
import os, sys, math, subprocess, pyfits, ConfigParser, numpy
from optparse import OptionParser
from StringIO import StringIO

# Check if parallel python module is present
try:
    import pp
except:
    print >> sys.stderr, 'Error: Parallel python module not found. Exiting.'
    sys.exit(-1)


# Load required IRAF packages
def loadPackages():
    pyraf.iraf.stsdas(_doprint = 0)
    pyraf.iraf.stsdas.analysis(_doprint = 0)
    pyraf.iraf.stsdas.analysis.restore(_doprint = 0)
    
    pyraf.iraf.noao(_doprint = 0)
    pyraf.iraf.noao.artdata(_doprint = 0)
    pyraf.iraf.noao.digiphot(_doprint = 0)
    pyraf.iraf.noao.digiphot.daophot(_doprint = 0)

    
# Generate x,y pixel positions for psf generation
def getPsfCoords(dimx, dimy, radius = 64):
    print >> sys.stdout, '\n Generating pixel positions for psf generation...'
    
    xy = []
    i = radius
    while i < dimx:
        j = radius
        while j < dimy:
            if i % 128 == 0 and j % 128 == 0:
                xy.append( (i, j) )
            j += radius
        i += radius    
    
    return xy    


# Generate image sections to be deconvolved. 256x256 sections are taken
def getChunks(dimx, dimy, radius = 64):
    psfxy = getPsfCoords(dimx, dimy, radius)
    
    print >> sys.stdout, '\n Generating image sections for deconvolution...'
    
    in_coords = []
    for value in psfxy:
        mix1, mix2 = value[0] - 127, value[0] + 128
        miy1, miy2 = value[1] - 127, value[1] + 128
            
        if mix1 < 1 and mix1 > -1e+8:
            mix1 = 1
        else:
            mix1 = mix1

        if miy1 < 1 and miy1 > -1e+8:
            miy1 = 1
        else:
            miy1 = miy1                
        
        if mix2 > dimx:
            mix2 = dimx
        else:
            mix2 = mix2
                
        if miy2 > dimy:
            miy2 = dimy
        else:
            miy2 = miy2                
        
        in_coords.append((mix1, mix2, miy1, miy2))
    
    return psfxy, in_coords


# Create input and output coordinates for pasting
def pasteCoords(imgxy, dimx, dimy, nsub):
    sub_coords = []
        
    for value in imgxy:
        sub_coords.append((nsub * (value[0] - 1) + 1, nsub * value[1], nsub * (value[2] - 1) + 1, nsub * value[3]))
    
    mem_coords, out_coords = [], []
    # Generate MEM output coordinates
    for value in sub_coords:
        # X1 and X2 values
        if value[0] != 1:
            x1 = value[0] + 64 * nsub
        else:
            x1 = value[0]
                
        if value[1] != (dimx * nsub):
            x2 = value[1] - 64 * nsub
        else:
            x2 = value[1]
                
        # Y1 and Y2 values
        if value[2] != 1:
            y1 = value[2] + 64 * nsub
        else:
            y1 = value[2]
                
        if value[3] != (dimy * nsub):
            y2 = value[3] - 64 * nsub
        else:
            y2 = value[3]
                
        out_coords.append((x1, x2, y1, y2))    
    
    # Generate MEM coordinates for input deconvolved sections            
    for value in out_coords:
        # X1 abd X2 coordinates
        if value[0] != 1 and value[0] > -1e+8:
            x1 = ( 64 * nsub ) + 1
        elif value[0] > -1e+8:
            x1 = 1
        else:
            x1 = -1e+16
                
        if value[1] != nsub * dimx and value[1] > -1e+8:
            x2 = 192 * nsub   
        elif value[1] > -1e+8:
            x2 = 160 * nsub
        else:
            x2 = -1e+16
                
        # Y1 and Y2 coordinates
        if value[2] != 1 and value[2] > -1e+8:
            y1 = ( 64 * nsub ) + 1
        elif value[2] > -1e+8:
            y1 = 1
        else:
            y1 = -1e+16
                
        if value[3] != nsub * dimx and value[3] > -1e+8:
            y2 = 192 * nsub   
        elif value[3] > -1e+8:
            y2 = 160 * nsub
        else:
            y2 = -1e+16            
                
        mem_coords.append((x1, x2, y1, y2))
    
    return mem_coords, out_coords           


# Function to copy FITS image or image section (uses pyfits to handle FITS files)
def imcopy(infile, outfile, dim = None):
    print >> sys.stdout, 'Copying ', infile, ' ----> ', outfile
    
    if len(outfile.split('[')) == 1:
        subprocess.call('cp ' + infile + '  ' + outfile, shell = True)
    else:
        if not dim:
            print >> sys.stderr, 'Error : for image section copying, dim parameter cannot be None. Exiting.'
            sys.exit(-1)
            
        header = pyfits.getheader(infile)
        output = numpy.zeros((dim, dim), dtype = numpy.float32)
        
        try:
            f1 = pyfits.open(infile)
        except:
            print >> sys.stderr, 'Error : Not able to open ', infile, '. Exiting.'
            sys.exit(-1)
    
        x1, x2 = int(outfile.split('[')[1].replace(']', '').split(',')[0].split(':')[0] ), int(outfile.split('[')[1].replace(']', '').split(',')[0].split(':')[1])
        y1, y2 = int(outfile.split('[')[1].replace(']', '').split(',')[1].split(':')[0] ), int(outfile.split('[')[1].replace(']', '').split(',')[1].split(':')[1])
        output[x1:x2, y1:y2] = f1[0].data

        outfile = outfile.split('[')[0]
        subprocess.call('rm -f ' + outfile, shell = True)
        pyfits.writeto(outfile, output, header = header)
        
    return outfile


# Generate psf image from psf lookup tables
def seepsf(psf, x, y, dim):
    seepsfimg_t = psf.replace('.fits', '_seepsf_' + str(x) + '_' + str(y) + 't.fits')
    subprocess.call('rm -f ' + seepsfimg_t, shell = True)
    
    print >> sys.stdout, '\n Generating psf image : ', seepsfimg_t
     
    pyraf.iraf.seepsf(psf, seepsfimg_t, dimension = dim, x = x, y = y)

    n = math.log(dim, 2)
    
    if math.ceil(n) != math.floor(n):
        psfdim = int(2**(math.floor(n) + 1))
        seepsfimg = psf.replace('.fits', '_seepsf_' + str(x) + '_' + str(y) + '.fits') 
        
        # Using imcopy method instead of IRAF imcopy function as IRAF imcopy task in failing in parallel mode
        if (psfdim - dim) % 2 == 0:
            cut = (psfdim - dim) / 2
            imcopy(seepsfimg_t, seepsfimg + '[' + str(cut) + ':' + str(psfdim - cut) + ',' + str(cut) + ':' + str(psfdim - cut) + ']', psfdim)
        else:
            cut1 = (psfdim - dim) / 2
            cut2 = (psfdim - dim) - cut1
            imcopy(seepsfimg_t, seepsfimg + '[' + str(cut2) + ':' + str(psfdim - cut1) + ',' + str(cut2) + ':' + str(psfdim - cut1) + ']', psfdim)
    else:
        seepsfimg = psf.replace('.fits', '_seepsf_' + str(x) + '_' + str(y) + '.fits') 
        subprocess.call('rm -f ' + seepsfimg , shell = True)
        subprocess.call('cp ' + seepsfimg_t + ' ' + seepsfimg, shell = True)
     
    # Cleanup - remove temprary seepsf files
    subprocess.call('rm -f ' + seepsfimg_t, shell = True)
            
    return seepsfimg


# Create output blank image
def createBlankImage(image, outfile = None, nsub = 1):
    if not outfile:
        if len(image.rsplit('[', 1)) > 1:
            outfile = image.rsplit('[', 1)[0].replace('.fits', '.' + image.rsplit('[', 1)[1].replace(']', '') + '_mem.fits')    
        else:
            outfile = image.rsplit('[', 1)[0].replace('.fits', '.mem.1.fits')
    
    subprocess.call('rm -f ' + outfile, shell = True)

    pyraf.iraf.blkrep(image, outfile, b1 = nsub, b2 = nsub)

    pyraf.iraf.imreplace(outfile, 0, upper = 'INDEF', lower = 'INDEF')

    return outfile


# Deconvolve image section using Maximum Entropy Method
def mem(image, psf, psfxy, imgxy, psfrad, nsub):
    # Load IRAF packages
    loadPackages()
    
    # Generate psf image for the image section
    dim = 2 * int(psfrad) * nsub + 1
    
    # Create psf image from psf tables
    psfimg = seepsf(psf, psfxy[0], psfxy[1], dim)

    outimg = image.rsplit('[', 1)[0].replace('.fits', '.mem_' + str(psfxy[0]) + '_' + str(psfxy[1]) + '.fits')
    subprocess.call('rm -f ' + outimg, shell = True)

    # Read MEM deconvolution parameters from configuration file
    parser = ConfigParser.SafeConfigParser()
    if not parser.read('deconvolve.cfg'):
        print >> sys.stderr, 'Error: Not able to open deconvolve.cfg configuraton file. Exiting.'
        sys.exit(-1)
        
    # Deconvolve the image section
    pyraf.iraf.mem(image + "[" + str(imgxy[0]) + ":" + str(imgxy[1]) + "," + str(imgxy[2]) + ":" + str(imgxy[3]) + "]", psf = psfimg, model = '', output = outimg, noise = parser.get('mem', 'noise'), adu = parser.get('mem', 'adu'), nsub = nsub, poisson = parser.get('mem', 'poisson'), tp = parser.get('mem', 'tp'), blksum = parser.get('mem', 'blksum'), guess = parser.get('mem', 'guess'), icf = parser.get('mem', 'icf'), hidden = parser.get('mem', 'hidden'), aim = parser.get('mem', 'aim'), maxiter = parser.get('mem', 'maxiter'), message = parser.get('mem', 'message'), m_update = parser.get('mem', 'm_update'), method = parser.get('mem', 'method'))
    
    # Clean up - remove seepsf file
    subprocess.call('rm -f ' + psfimg, shell = True)

    return outimg


# Combine overlapping deconvolved image secitons to generate final deconvolved image
def paste(image, mem_imgs, imgxy, dimx, dimy, nsub, outfile):
    memimg = createBlankImage(image, outfile, nsub)
    
    # Determine coordinates for pasting the image
    mem_coords, out_coords = pasteCoords(imgxy, dimx, dimy, nsub)
 
    # Paste image sections to create the final deconvolved image
    for i in range(len(mem_coords)):
       pyraf.iraf.imcopy(mem_imgs[i] + '[' + str(mem_coords[i][0]) + ':' + str(mem_coords[i][1]) + ',' + str(mem_coords[i][2]) + ':' + str(mem_coords[i][3]) + ']', memimg + '[' + str(out_coords[i][0]) + ':' + str(out_coords[i][1]) + ',' + str(out_coords[i][2]) + ':' + str(out_coords[i][3]) + ']')
    
    # Clean up - remove the MEM deconvolved image sections
    for i in range(len(mem_coords)):
        subprocess.call('rm -f ' + mem_imgs[i], shell = True)
                
    return memimg
    

# Parallel python worker method
# =============================
def worker(indata):
    cnt, image, psf, psfxy, imgxy, psfrad, nsub = indata
    
    result = mem(image, psf, psfxy, imgxy, psfrad, nsub)
    
    return (cnt, result)


# Main function of deconvolution routine
# ======================================
def deconvolve(image, psf, nsub, outfile = None, ncpus = None):
    # Determine image dimensions
    pyraf.iraf.imgets(image, 'i_naxis1')
    dimx = int(pyraf.iraf.imgets.value)
    
    pyraf.iraf.imgets(image, 'i_naxis2')
    dimy = int(pyraf.iraf.imgets.value)

    # Determine PSF radius
    pyraf.iraf.imgets(psf, 'psfrad')
    psfrad = float(pyraf.iraf.imgets.value)
    
    # Calculate image section coordinates
    psfxy, imgxy = getChunks(dimx, dimy)

    # Get number of processors (cores) on the machine. In case of Intel processors with
    # hyperthreading, total number of processors will be equal to number of cores * number 
    # of threads/core. Default is maximum number of cores available.
    ppservers = ()
    if ncpus:
        # Creates jobserver with ncpus workers
        job_server = pp.Server(int(ncpus), ppservers=ppservers)
    else:
        # Creates jobserver with automatically detected number of workers
        job_server = pp.Server(ppservers = ppservers)

    # Generate dataset chunks to be processed
    jobs = []
    for i in range(len(psfxy)):
        chunks = (i, image, psf, psfxy[i], imgxy[i], psfrad, nsub) 
        jobs.append(job_server.submit(worker, (chunks,), (mem,seepsf,loadPackages,imcopy,), ("sys","ConfigParser","pyraf","subprocess","pyfits","numpy","math",)))

    # wait for all the jobs to finish
    job_server.wait()

    # Append the results
    results = []
    for job in jobs:
        results.append(job())

    # Sort the result list to have right series for pasting
    results = sorted(results, key = lambda result: result[0])
    
    # Generate final deconvolved image by pasting image sections
    mem_imgs = []
    for value in results:
       mem_imgs.append(value[1])
        
    deconv_img = paste(image, mem_imgs, imgxy, dimx, dimy, nsub, outfile)
    
    print >> sys.stdout, '\n Final Deconvolved Image : ', deconv_img, '\n'



# Input validation funtion, calls main deconvolution routine
# ==========================================================
def main(image, psf, nsub, outfile = None, ncpus = None):
    # Check if the image exists
    if not os.path.exists(image.rsplit( '[', 1 )[0]):
        print >> sys.stderr, 'Error: Image ', image, ' does not exist. Exiting.'
        sys.exit(-1)

    # If multi-extension file, check if the user entered the extension or not
    if len( pyfits.open( image.rsplit( '[', 1 )[0] ) ) > 1 and len( image.rsplit( '[', 1 ) ) == 1:
        print >> sys.stdout, 'Error : Multi-extension FITS image. Please provide image extension. Exiting.'
        sys.exit(-1)

    # Check if psf file exists
    if not os.path.exists( psf ):
        print >> sys.stderr, 'Error: PSF file ', psf, ' does not exist. Exiting.'
        sys.exit(-1)	

    # Execute the method
    deconvolve(image, psf, int( nsub ), outfile, ncpus)



# Entry point for deconvolution utility
# =====================================
if __name__ == '__main__':
    usage = "Usage: python %prog [options] image psf nsub"
    description = "Description. Utility to deconvolve image with spatially varying psf in multiprocessing mode.\nMaximum entropy method is used with subsampling (nsub should be greater than 1)."
    parser = OptionParser(usage = usage, version = "%prog 1.0", description = description)
    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose", default = False,
                      help = "print result messages to stdout"
                      )
    parser.add_option("-q", "--quiet",
                    action="store_false", dest="verbose", default = True,
                    help = "don't print result messages to stdout"
                    )    
    parser.add_option("-o", "--outfile", dest = "filename",
                    action='store', metavar="FILE", help = "output file name"
                    )
    parser.add_option("-n", "--ncpus", dest = "ncpus", metavar="NCPUS",
                    action="store", help = "number of cpus (cores) for processing"
                    )
    (options, args) = parser.parse_args()
    
    # Check for number of input arguments
    if len(args) != 3:
        parser.error("Incorrect number of arguments. Check help for more details.")

    if args[2] <= 1:
        parser.error("Subsampling should be greater than 1. Check help for more details.")
 
        
    print >> sys.stdout, '\n Starting processing...'
    
    # Check verbosity
    if not options.verbose:
        output = StringIO()
        old_stdout = sys.stdout
        sys.stdout = output

    # Check if pyraf module is installed
    try:
        import pyraf
    except:
        print >> sys.stderr, 'Error: Python module pyraf not found. Exiting.'
        sys.exit(-1)
    
    main(args[0], args[1], args[2], options.filename, options.ncpus)
    
    # Reset verbosity
    if not options.verbose:
        sys.stdout = old_stdout
    
    print >> sys.stdout, '\n Process completed successfully.'