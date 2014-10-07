#!/usr/bin/env python
""" Analyse GALFACTS Stokes image cubes of calibrator sources for a given field.
    Russ Taylor, 22 January 2014
    Modified August 2014 to include corrections for absolute K/Jy in each beam.
"""

import pyfits
from math import *
from scipy import *
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import sys, os, time, shutil

from calibdefs import *            # classes and functions for calib.py

#----------------------------------------------------------------------
# initialise constants and arrays
rad2deg = 180.0/3.1415928
deg2rad = 1.0/rad2deg
latitude = 18.3442               # latitude of Arecibo Observatory
c = 3.e8                         # speed of light (m/s)
k = 1.38e-23                     # Boltzmann constant (J/K)
bmaj = 4.0                       # beam size in arc minutes
solid_angle = (pi/(4.0*log(2)))*((bmaj/60.0)*(pi/180.0))**2   # beam solid angle in steradians
StoT = (c**2/(2.0*k)) * ((1e-06)**2) * (1e-26) / solid_angle
#print StoT

startchan = 263
HIchans = (2448,2528)
smoothwidth = 100
smooth = True
debug = False
plotbeamcorrs = False     # set true if I want one plot per beam for the corrections
poltol = 0.01             # fractional polarization upper limit for assumed unpolarized

alldecs = []       # array for declinations of all sources over all fields
alldecs6= []       # array for declinations of all sources for beam 6 over all fields
allelev = []       # array for elevations of all sources over all fields
allKpJy = []
allrelgain = []

allKpJ0 = []       # array for all beam 0 gains of all sources over all fields
allKpJ1 = []       # array for all beam 1 gains of all sources over all fields
allKpJ2 = []       # array for all beam 2 gains of all sources over all fields
allKpJ3 = []       # array for all beam 3 gains of all sources over all fields
allKpJ4 = []       # array for all beam 4 gains of all sources over all fields
allKpJ5 = []       # array for all beam 5 gains of all sources over all fields
allKpJ6 = []       # array for all beam 6 gains of all sources over all fields

relgain1 = []       # array for all beam 1 gains of all sources over all fields
relgain2 = []       # array for all beam 2 gains of all sources over all fields
relgain3 = []       # array for all beam 3 gains of all sources over all fields
relgain4 = []       # array for all beam 4 gains of all sources over all fields
relgain5 = []       # array for all beam 5 gains of all sources over all fields
relgain6 = []       # array for all beam 6 gains of all sources over all fields

sr = [50,100,50,100]      # region containing the source [x1,x2,y1,y2]
soff = [20,120,110,120]   # region for calculating the off source brightness [x1,x2,y1,y2]

#fields = ['S1','S2','S3','S4','N1','N2','N3','N4']
#fields = ['S2','S3','S4','N2','N3','N4']
fields = ['N4']
#--------------------------------------------------------------------------
for field in fields:
    beamnumbers = ['0','1','2','3','4','5','6'] 
    if field in ['S2','S3','S4']:
        beamnumbers = ['0','1','2','3','4','5'] 
    nbeams = len(beamnumbers)
    freq = []                         # array of frequencies for each channel
    solutions = []                    # list for beam solutions for each source
    KpJy = []                         # array for K per Jy per beam per source
    relgains = []                     # array for gains relative to beam 0 per source
    
    sources =[]
    indir = field+'/images/' 
    sourcefilename = field+'/'+field+'.sources.txt'
    sourcefile = open(sourcefilename,'r')
    sourcelines = sourcefile.readlines()
    sourcefile.close()
#    print "\n Read %d lines from %s" % (len(sourcelines), sourcefilename)
    for j in range(len(sourcelines)):
        name = sourcelines[j].split()[0]
        if name[0] != '#':
            RA = float(sourcelines[j].split()[1])
            DEC = float(sourcelines[j].split()[2])
            Flux = float(sourcelines[j].split()[3])
            Index = float(sourcelines[j].split()[4])
            Pflux = float(sourcelines[j].split()[5])
            PPA = float(sourcelines[j].split()[6])
            A1 = float(sourcelines[j].split()[7])
            A2 = float(sourcelines[j].split()[8])
            A3 = float(sourcelines[j].split()[9])
            sources.append( source(name,RA,DEC,Flux,Index,Pflux,PPA,(A1,A2,A3)) )

    print "\nUsed %2d of %2d sources in %s" % (len(sources), len(sourcelines), sourcefilename)
    relativegain = []    # array for mean relative beam gains over all sources
    relativegainrms = []  # array for rms of relative beam gains over all sources
    decs = []       # array for declinations of calibration sources 
    for s in sources:
        decs.append(s.dec)

    outdir = field+'/'+field+'.calib.true.out/'
    try:
        shutil.rmtree(outdir)
        print " deleting existing output directory: %s" % outdir
    except: print " output directory does not exist"
    print " creating new output directory: %s\n" % outdir
    os.mkdir(outdir)

# determine frequencies for each channel

    Qcorrfile = open(outdir+field+'.Qcorr.txt','w')
    Ucorrfile = open(outdir+field+'.Ucorr.txt','w')
    Vcorrfile = open(outdir+field+'.Vcorr.txt','w')
    Icorrfile = open(outdir+field+'.Icorr.txt','w')
    Beamgainfile = open(outdir+field+'.Beamgain.txt','w')
    outfile = open(outdir+field+'.calib.out','w')
 
    Ifile = indir+sources[0].name+'_BEAM0_0263_4022_I.fits'
    print "\n---------------"

    hdulist = pyfits.open(Ifile)

    channels        = hdulist[0].header['NAXIS3']
    ReferenceFreq   = hdulist[0].header['CRVAL3']
    channelwidth    = hdulist[0].header['CDELT3']
    ReferenceChannel = hdulist[0].header['CRPIX3']

    for i in range(channels):
        frequency = (ReferenceFreq + (i+1 - ReferenceChannel)*channelwidth)/1e06
        freq.append( frequency )
    freq = array(freq)
    meanfreq = freq.mean()
    bandwidth = fabs(freq[len(freq)-1] - freq[0])

    print "Processing %s calibrators" % field
    print "Number of beams = %2d" % nbeams
    print "Number of channels = %d" % channels
    print "Mean Frequency = %8.2f MHz" % (meanfreq)
    print "Bandwidth = %8.2f MHz" % bandwidth

    outfile.write("Processing %s calibrators\n" % field)
    outfile.write("Number of channels = %d\n" % channels)
    outfile.write("Mean Frequency = %8.2f MHz\n" % meanfreq)
    outfile.write("Bandwidth = %8.2f MHz\n " % bandwidth)

    nvsscorr = (1400/meanfreq)**(0.7)
    hdulist.close()

#--------------------------------------------------------
# solve for beam solutions for each source
    for s in sources:
        print "\nProcessing source %s" % s.name 
        outline = "Processing source %s \n" % s.name 
        outfile.write(outline)
        outfile.write(outline)
        srcflux = []                 # real flux density of source in each channel
        srctemp = []                 # brightness temperature of the source in each channel
        beams = []                   # array of beam properties for this source 
        gains = []                   # array of gains in K per Jy for this source

 # get source total intensity spectrum   
 
        for i in range(channels):
            srcflux.append( fluxdensity(s.a,freq[i]) )
            srctemp.append( fluxdensity(s.a,freq[i])/(freq[i]**2))
        srcflux = array(srcflux)
        srctemp = array(srctemp)*StoT
# normalize the Salter solution so that it gives the NVSS flux density at 1400 MHz.        
        S1400 = fluxdensity(s.a,1400)
        fluxcorrection = s.flux/S1400
        sourceflux = srcflux.mean()*fluxcorrection/1000.0       

        meansrcfracpol = s.P/s.flux
        KperJy = srctemp.mean()/srcflux.mean()
        srcfluxnorm = srcflux/srcflux.mean()
        srctempnorm = srctemp/srctemp.mean()
       
        print "   Mean flux density over the band = %6.3f Jy, %6.2f K" %(srcflux.mean(), srctemp.mean())
        print "   Corr flux density over the band = %6.3f Jy." %( sourceflux)
        print "   NVSS flux density at mid band   = %6.3f Jy." % (nvsscorr*s.flux/1000.0)
        print "   NVSS Percent polarization = %7.2f" % (100.0*meansrcfracpol)
        print "   NVSS polarization position angle = %6.2f" % s.PA
    
        outline = "   Mean flux density over the band = %6.3f Jy, %6.2f K \n" %(srcflux.mean(), srctemp.mean())
        outfile.write(outline)
        outline = "   Corr flux density over the band = %6.3f Jy." % sourceflux
        outfile.write(outline)
        outline = "   NVSS Percent polarization = %7.2f\n" % (100.0*meansrcfracpol)
        outfile.write(outline)
        outline = "   NVSS polarization position angle = %6.2f \n" % s.PA
        outfile.write(outline)

    
# ---------- plot 
#    plt.plot(freq,srcflux)
#    plt.xlabel('frequency (MHz)')
#    plt.ylabel('Flux Density (Jy)')
#    plt.title(s.name)
#    plt.show()
 
# now loop through beams to derive corrections for each beam     
        for b in beamnumbers: 
            peakI = []; peakQ = []; peakU = []; peakV = []
            offI  = []; offQ  = []; offU  = []; offV  = []
            polPA = []; PAcorr = []
            gaussI = []; gaussQ = []; gaussU= []; gaussV =[]
            bxI = []; byI = []; bpaI = []; gaussx0= []; gaussy0 = []
 
 #  load up the array for the full resolution spectral correction for this beam 
            
            Ifile = indir+s.name+'_BEAM'+b+'_0263_4022_I.fits'
            Qfile = indir+s.name+'_BEAM'+b+'_0263_4022_Q.fits'
            Ufile = indir+s.name+'_BEAM'+b+'_0263_4022_U.fits'
            Vfile = indir+s.name+'_BEAM'+b+'_0263_4022_V.fits'
        
            hdulistI = pyfits.open(Ifile)
            hdulistQ = pyfits.open(Qfile)
            hdulistU = pyfits.open(Ufile)
            hdulistV = pyfits.open(Vfile)

            Idata = hdulistI[0].data
            Qdata = hdulistQ[0].data
            Udata = hdulistU[0].data
            Vdata = hdulistV[0].data
            
            hdulistI.close()
            hdulistQ.close()
            hdulistU.close()
            hdulistV.close()
    
            nx = Idata.shape[2]
            ny = Idata.shape[1]
            nchan = Idata.shape[0]
            chan = arange(0,nchan,1)

            x = linspace(0,nx,ny)    # uses for experimenta gaussian fit version
            y = linspace(0,nx,ny)
            x,y = meshgrid(x,y)


 # ----------------------------------
# now analyse image cubes
    
# find pixel coordinates of the data maximum in I in channel 20           
            maxindex = Idata[20,:,:].argmax()
            ay = Idata[20,sr[2]:sr[3],sr[0]:sr[1]].max(axis=0)
            ax = Idata[20,sr[2]:sr[3],sr[0]:sr[1]].max(axis=1)
            xmax = sr[0] + ay.argmax()
            ymax = sr[2] + ax.argmax()         
            if debug:
                print Idata[20,:,:].shape
                print ax
                print ay
                print ax.argmax()
                print ay.argmax()
                print "\n Beam %s, (xmax,ymax) = %d, %d" % (b, xmax,ymax)
                print " value = %f" % Idata[20,ymax,xmax]
                print Idata[i,:,:].shape

#  this code is for fitting gaussian to the beam
#        for i in range (nchan):
#            initial_guess = (Idata[i,ymax,xmax],xmax,ymax,3.0,3.0,0.0,0.0)
#            popt,pcov = curve_fit(twoDGaussian,(x,y),Idata[i,:,:].ravel(),p0=initial_guess)
#            print popt
#            gaussI.append(popt[0])
#            bxI.append(popt[1])
#            byI.append(popt[2])
#            bpaI.append(popt[3])

# -- now find I, Q, U and V data values a position of I peak and remove background to get
# spectra of I, Q, U and V at the position of the I peak
            for i in range(nchan):
                peakI.append(Idata[i,ymax,xmax])
                peakQ.append(Qdata[i,ymax,xmax])
                peakU.append(Udata[i,ymax,xmax])    
                peakV.append(Vdata[i,ymax,xmax]) 
                offI.append( Idata[i,soff[2]:soff[3],soff[0]:soff[1]].mean() )
                offQ.append( Qdata[i,soff[2]:soff[3],soff[0]:soff[1]].mean() )
                offU.append( Udata[i,soff[2]:soff[3],soff[0]:soff[1]].mean() )
                offV.append( Vdata[i,soff[2]:soff[3],soff[0]:soff[1]].mean() )

            peakI = removeRFI(array(peakI))
            peakI = removeRFI(array(peakI))
            peakQ = removeRFI(array(peakQ))
            peakU = removeRFI(array(peakU))
            peakV = removeRFI(array(peakV))
            offI = removeRFI(array(offI))
            offQ = removeRFI(array(offQ))
            offU = removeRFI(array(offU))
            offV = removeRFI(array(offV))           
            Isrc = removeRFI(peakI - offI)
            Qsrc = removeRFI(peakQ - offQ)
            Usrc = removeRFI(peakU - offU)
            Vsrc = removeRFI(peakV - offV)
        
# fix channels in HI line region
            value = (Isrc[HIchans[0]-1] + Isrc[HIchans[1]+1])/2.0        
            for i in range(HIchans[0],HIchans[1]):
               Isrc[i]=value
            
            meanI = Isrc.mean()          # mean value of Stokes I spectrum
            normI = Isrc/meanI           # normalize I spectrum to mean of 1.0   
            Icorr = srctempnorm/normI    # correction to get true source temperature spectrum
            Itrue = Isrc*Icorr
            r = rms(Icorr[200:500])      # rms of correction spectrum from chan 200 to 500
        
            if(smooth):
                Icorr = boxcar(smoothwidth,Icorr)
                Isrc  = boxcar(smoothwidth,Isrc)
                Qsrc  = boxcar(smoothwidth,Qsrc)
                Usrc  = boxcar(smoothwidth,Usrc)
                Vsrc  = boxcar(smoothwidth,Vsrc)      
            fracQ = Qsrc/Itrue
            fracU = Usrc/Itrue
            fracV = Vsrc/Itrue
    
            fracP = sqrt(fracQ*fracQ + fracU*fracU)   # spectrum of fractional polarizaiton
            Pcorr = meansrcfracpol/fracP              # correction for observed to true P       
            for i in range(len(fracU)):
                polPA.append(rad2deg*0.5*atan2(fracU[i],fracQ[i]))
            for i in range(len(polPA)-1):
                if ((polPA[i+1] - polPA[i]) < -100 ):
                    polPA[i+1] = polPA[i+1] + 180.0
            for i in range(len(polPA)):
                PAcorr.append(s.PA-polPA[i])
            PAcorr = array(PAcorr)
        
#  derive mean values over the band
            meanFracQ = Qsrc.mean()/meanI
            meanFracU = Usrc.mean()/meanI
            meanFracV = Vsrc.mean()/meanI
            meanFracP = sqrt(meanFracQ*meanFracQ + meanFracU*meanFracU)
            meanPA = rad2deg* 0.5 * atan2(meanFracU,meanFracQ)
            meanPcorr = meansrcfracpol/meanFracP
            meanPAcorr =  (s.PA - meanPA)
    
#---------------------------------   
#            plt.plot(chan,Itrue)
#            plt.plot(chan,Isrc)
#            plt.title('Beam: ' + b)
#            plt.show()

# print and write out results       
            print "Beam %s: I =%7.3f, Q = %7.3f, U = %7.3f, V = %7.3f, K/Jy = %7.2f" % \
                 (b, Isrc.mean(),Qsrc.mean(),Usrc.mean(),Vsrc.mean(), meanI/(sourceflux))
            print "                    Q/I=%7.3f, U/I=%7.3f, V/I=%7.3f, P=%6.2f PA=%6.1f, DPA = %5.1f, rms = %6.3f" % \
                 (meanFracQ, meanFracU, meanFracV, 100.0*meanFracP, meanPA, meanPAcorr, r)     
             
            outline =  "Beam %s: I =%7.3f, Q = %7.3f, U = %7.3f, V = %7.3f, K/Jy = %7.2f\n" % \
                 (b, Isrc.mean(),Qsrc.mean(),Usrc.mean(),Vsrc.mean(), meanI/(sourceflux))  
            outfile.write(outline)   
            outline =  "                    Q/I=%7.3f, U/I=%7.3f, V/I=%7.3f, P=%6.2f PA=%6.1f,DPA = %5.1f,  rms = %6.3f\n" % \
                 (meanFracQ, meanFracU, meanFracV, 100.0*meanFracP, meanPA, meanPAcorr, r)       
            outfile.write(outline)  

            gains.append(meanI/(sourceflux))  
            beams.append( beam(b, Isrc, fracQ, fracU, fracV, polPA, Icorr, r, Pcorr, PAcorr) )
     
        KpJy.append(gains)     
        solutions.append(beams) # add the beam properties for this source to the solution list
        outfile.write("\n")
        relgains.append(gains/gains[0])
        alldecs.append(s.dec)
        allKpJ0.append(gains[0])
        allKpJ1.append(gains[1])
        allKpJ2.append(gains[2])
        allKpJ3.append(gains[3])
        allKpJ4.append(gains[4])
        allKpJ5.append(gains[5])
        relgain1.append(gains[1]/gains[0])
        relgain2.append(gains[2]/gains[0])
        relgain3.append(gains[3]/gains[0])
        relgain4.append(gains[4]/gains[0])
        relgain5.append(gains[5]/gains[0])
        if (nbeams==7):
            allKpJ6.append(gains[6])
            relgain6.append(gains[6]/gains[0])
            alldecs6.append(s.dec)
#        print relgain1
# now average over all sources to create an average correction spectrum for each beam

#print KpJy

#------------K per Jy and gains----------------------------
    beamgain = []
    relbeamgain = []
    for j in range(len(beamnumbers)):
        temp1 = []; temp2=[]
        for s in KpJy:
            temp1.append(s[j])
        for s in relgains:
            temp2.append(s[j])
        beamgain.append(temp1)
        relbeamgain.append(temp2)
    
    beamgain=array(beamgain)
    relbeamgain=array(relbeamgain)
    beam0KpJy = beamgain[0].mean()
    beam0KpJyerr = beamgain[0].std()


    print "\n Average beam 0 gain = %8.3f +/- %7.3f" % (beam0KpJy, beam0KpJyerr)
    outline = "\n Average beam 0 gain = %8.3f +/- %7.3f\n" % (beam0KpJy, beam0KpJyerr)
    outfile.write(outline)

    for j in range(len(beamnumbers)):
        print "%s %8.4f %7.4f" % (beamnumbers[j], relbeamgain[j].mean(), relbeamgain[j].std())
        outline= "%s %8.4f %7.4f \n" % (beamnumbers[j], relbeamgain[j].mean(), relbeamgain[j].std())
        outfile.write(outline)
        Beamgainfile.write(outline)
        relativegain.append( relbeamgain[j].mean() )
        relativegainrms.append( relbeamgain[j].std() )
    Beamgainfile.close()
        
    plt.figure(figsize = (10,8))
    for j in range(len(beamnumbers)):
       plt.plot(decs, beamgain[j],label = 'Beam: '+beamnumbers[j])
    plt.ylim([5,15])  
    plt.xlim([-2.0,18.0])
    plt.legend(loc='best')
    plt.xlabel('Declination')
    plt.ylabel('Beam Gain (K per Jy)')
    plt.savefig(outdir+'Beamgains')
#   plt.show()
    plt.close()

    plt.figure(figsize = (10,8))
    for j in range(len(beamnumbers)):
       plt.plot(decs, relbeamgain[j],label = 'Beam: '+beamnumbers[j])
    plt.ylim([0.6,1.1])
    plt.xlim([-2.0, 18.0])
    plt.legend(loc='best')
    plt.xlabel('Declination')
    plt.ylabel('Relative beam gains')
    plt.savefig(outdir+'RelBeamgains')
#   plt.show()
    plt.close()


#--------I-----------
    beamcorr = []                   # array of beam solutions 
    for j in range(len(beamnumbers)):
        sumweights= 0
        temp = zeros( ( len(solutions[0][0].Icorr) ), float)  
        for s in solutions:
            if (s[j].rms < 0.1):
                w = 1.0/(s[j].rms*s[j].rms)     # weight for the solution for this source
                sumweights = sumweights + w
                for i in range(len(s[j].Icorr)):
                    temp[i] = temp[i] + w * s[j].Icorr[i]   
        beamcorr.append(temp/sumweights)
    
    plt.figure(figsize=(10,8))
    plt.ylim([0.5,1.5])
    for j in range(len(beamnumbers)):
        plt.plot(chan,beamcorr[j],label = 'Beam: '+beamnumbers[j])  
    plt.xlabel('Channel') 
    plt.ylabel('Spectral correction')   
    plt.legend(loc='best')
    plt.title('Noise weighted average Stokes I spectral correction of '+field+' calibrators') 
    plt.savefig(outdir+'AveCorrections')
    plt.close()
#   plt.show() 

# ------relative gains in I ------------------------
    plt.figure(figsize=(10,8))
    plt.ylim([0.5,1.5])
    for j in range(len(beamnumbers)):
        plt.plot(chan,beamcorr[j]/relativegain[j],label = 'Beam: '+beamnumbers[j])  
    plt.xlabel('Channel') 
    plt.ylabel('Spectral correction')   
    plt.legend(loc='best')
    plt.title('Noise weighted average Stokes I spectral correction of '+field+' calibrators') 
    plt.savefig(outdir+'AveRelativeCorrections')
    plt.close()
#   plt.show() 
    
#--------V-----------    
    Vcorr = []                   # array of beam solutions 
    for j in range(len(beamnumbers)):
        sumweights= 0
        temp = zeros( ( len(solutions[0][0].V) ), float)  
        for s in solutions:
            if (s[j].rms < 0.1):
                w = 1.0/(s[j].rms*s[j].rms)     # weight for the solution for this source
                sumweights = sumweights + w
                for i in range(len(s[j].V)):
                    temp[i] = temp[i] + w * s[j].V[i]   
        Vcorr.append(temp/sumweights)
    
    plt.figure(figsize=(10,8))
    plt.ylim([-0.2,0.2])
    for j in range(len(beamnumbers)):
        plt.plot(chan,Vcorr[j],label = 'Beam: '+beamnumbers[j])  
    plt.xlabel('Channel') 
    plt.ylabel('V Leakage')   
    plt.legend(loc='best')
    plt.title('Noise weighted average of V leakage of '+field+' calibrators') 
    plt.savefig(outdir+'AveVCorrections')
    plt.close()

#-----------Q-----------------
    Qcorr = []                   # array of beam solutions 
    for j in range(len(beamnumbers)):
        sumweights= 0
        temp = zeros( ( len(solutions[0][0].Q) ), float) 
        count = 0 
        for s in solutions:
            if (s[j].rms < 0.1):
                sourcepol = sources[count].P/sources[count].flux
                if(sourcepol < poltol):
                    print "%s has polarization %5.3f. Using for Q leakage" % (sources[count].name, sourcepol)
                    w = 1.0/(s[j].rms*s[j].rms)     # weight for the solution for this source
                    sumweights = sumweights + w
                    for i in range(len(s[j].Q)):
                        temp[i] = temp[i] + w * s[j].Q[i] 
            count = count + 1  
        Qcorr.append(temp/sumweights)
    
    plt.figure(figsize=(10,8))
    plt.ylim([-0.2,0.2])
    for j in range(len(beamnumbers)):
        plt.plot(chan,Qcorr[j],label = 'Beam: '+beamnumbers[j])  
    plt.xlabel('Channel') 
    plt.ylabel('Q Leakage')   
    plt.legend(loc='best')
    titlestring = 'Noise weighted average of Q leakage of '+field+' calibrators'
    plt.title(titlestring) 
    plt.savefig(outdir+'AveQCorrections')
    plt.close()

#---------U---------------------
    Ucorr = []                   # array of beam solutions 
    for j in range(len(beamnumbers)):
        sumweights= 0
        temp = zeros( ( len(solutions[0][0].U) ), float)
        count = 0  
        for s in solutions:
            if (s[j].rms < 0.1):
                sourcepol = sources[count].P/sources[count].flux
                if(sourcepol < poltol):
                    print "%s polarization = %5.3f. Using for U leakage" % (sources[count].name, sourcepol)
                    w = 1.0/(s[j].rms*s[j].rms)     # weight for the solution for this source
                    sumweights = sumweights + w
                    for i in range(len(s[j].U)):
                        temp[i] = temp[i] + w * s[j].U[i]  
            count = count + 1 
        Ucorr.append(temp/sumweights)
    
    plt.figure(figsize=(10,8))
    plt.ylim([-0.2,0.2])
    for j in range(len(beamnumbers)):
        plt.plot(chan,Ucorr[j],label = 'Beam: '+beamnumbers[j])  
    plt.xlabel('Channel') 
    plt.ylabel('U Leakage')   
    plt.legend(loc='best')
    plt.title('Noise weighted average of U leakage of '+field+' calibrators')
    plt.savefig(outdir+'AveUCorrections')
    plt.close()
#   plt.show() 


#---------------------------------------  
# plot results per beam for all sources

    for j in range(len(beamnumbers)):
        plt.figure(figsize=(10,8))
        plt.ylim([0.5,1.5])
        for i in range(len(solutions)):
            plt.plot(chan,solutions[i][j].Icorr,label = sources[i].name)  
        plt.xlabel('Channel') 
        plt.ylabel('Bandpass Correction')   
        plt.legend(loc='best')
        plt.title(field+' calibrators: `Beam ' + beamnumbers[j])      
        plt.savefig(outdir+'Icorr_Beam'+beamnumbers[j])
        plt.close()

    for j in range(len(beamnumbers)):
        plt.figure(figsize=(10,8))
        plt.ylim([-0.2,0.2])
        for i in range(len(solutions)):
         plt.plot(chan,solutions[i][j].V,label = sources[i].name)  
        plt.xlabel('Channel') 
        plt.ylabel(' V fractional polarization')
        plt.legend(loc='best')
        plt.title(field+' calibrators: Beam ' + beamnumbers[j])      
        plt.savefig(outdir+'FracV_Beam'+beamnumbers[j])
        plt.close()
    
    for j in range(len(beamnumbers)):
        plt.figure(figsize=(10,8))
        plt.ylim([-0.2,0.2])
        for i in range(len(solutions)):
             plt.plot(chan,solutions[i][j].Q,label = sources[i].name)  
        plt.xlabel('Channel') 
        plt.ylabel(' Q fractional polarization')
        plt.legend(loc='best')
        plt.title(field+' calibrators: Beam ' + beamnumbers[j])      
        plt.savefig(outdir+'FracQ_Beam'+beamnumbers[j])
        plt.close()
    
    for j in range(len(beamnumbers)):
        plt.figure(figsize=(10,8))
        plt.ylim([-0.2,0.2])
        for i in range(len(solutions)):
             plt.plot(chan,solutions[i][j].U,label = sources[i].name)  
        plt.xlabel('Channel') 
        plt.ylabel(' U fractional polarization')
        plt.legend(loc='best')
        plt.title(field+' calibrators: Beam ' + beamnumbers[j])      
        plt.savefig(outdir+'FracU_Beam'+beamnumbers[j]) 
        plt.close()   

 
#--------------------------------------------------------
# create correction files with full 4096 channel range and plot for each beam, Stokes.
    Ispeccorr = []  
    Irelspeccorr = []  
    Qspeccorr = []
    Uspeccorr = []
    Vspeccorr = []
    fullres = range(4096)
    for j in range(len(beamnumbers)):
        Itemp= []; Qtemp = []; Utemp=[]; Vtemp= []; Ireltemp= []
        for i in range(4096):
            Itemp.append(1.0)
            Ireltemp.append(1.0)
            Qtemp.append(0.0)
            Utemp.append(0.0)
            Vtemp.append(0.0)
        for i in range(len(beamcorr[j])):
            Itemp[i+startchan] = beamcorr[j][i]
            Ireltemp[i+startchan] = beamcorr[j][i]/relativegain[j]
            Qtemp[i+startchan] = Qcorr[j][i]
            Utemp[i+startchan] = Ucorr[j][i]
            Vtemp[i+startchan] = Vcorr[j][i]
        Ispeccorr.append(Itemp)    
        Irelspeccorr.append(Ireltemp) 
        Qspeccorr.append(Qtemp)           
        Uspeccorr.append(Utemp)   
        Vspeccorr.append(Vtemp)      
        if(plotbeamcorrs):
            plt.figure(figsize=(8,6))
            plt.ylim([0.5,1.5])
            plt.title(outdir+'Icorr: Beam'+beamnumbers[j])
            plt.plot(fullres,Itemp)
            plt.savefig(outdir+'Beam'+beamnumbers[j]+'Icorr')
            plt.close()
        
            plt.figure(figsize=(8,6))
            plt.ylim([0.5,1.5])
            plt.title(outdir+'Irelcorr: Beam'+beamnumbers[j])
            plt.plot(fullres,Ireltemp)
            plt.savefig(outdir+'Beam'+beamnumbers[j]+'Irelcorr')
            plt.close()
    
            plt.figure(figsize=(8,6))
            plt.ylim([-0.2,0.2])
            plt.title(outdir+'Qcorr: Beam'+beamnumbers[j])
            plt.plot(fullres,Qtemp)
            plt.savefig(outdir+'Beam'+beamnumbers[j]+'Qcorr')
            plt.close()

            plt.figure(figsize=(8,6))
            plt.ylim([-0.2,0.2])
            plt.title(outdir+'Ucorr: Beam'+beamnumbers[j])
            plt.plot(fullres,Utemp)
            plt.savefig(outdir+'Beam'+beamnumbers[j]+'Ucorr')
            plt.close()
 
            plt.figure(figsize=(8,6))
            plt.ylim([-0.2,0.2])
            plt.title(outdir+'Vcorr: Beam'+beamnumbers[j])
            plt.plot(fullres,Vtemp)
            plt.savefig(outdir+'Beam'+beamnumbers[j]+'Vcorr')
            plt.close()
       
#------------------------------
# write out the correction files with 4096 channels.  
    if (nbeams == 6):
        for i in range(4096):
            outline = " %3d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f  1.0\n" \
             % (i, Ispeccorr[0][i], Ispeccorr[1][i], Ispeccorr[2][i], \
             Ispeccorr[3][i], Ispeccorr[4][i], Ispeccorr[5][i])   
            Icorrfile.write(outline)
        for i in range(4096):
            outline = " %3d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f  0.0\n" \
           % (i, Qspeccorr[0][i], Qspeccorr[1][i], Qspeccorr[2][i], \
            Qspeccorr[3][i], Qspeccorr[4][i], Qspeccorr[5][i])   
            Qcorrfile.write(outline)
        for i in range(4096):
            outline = " %3d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f  0.0 \n" \
            % (i, Uspeccorr[0][i], Uspeccorr[1][i], Uspeccorr[2][i], \
            Uspeccorr[3][i], Uspeccorr[4][i], Uspeccorr[5][i])#, Uspeccorr[6][i])   
            Ucorrfile.write(outline)
        for i in range(4096):
            outline = " %3d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f  0.0\n" \
           % (i, Vspeccorr[0][i], Vspeccorr[1][i], Vspeccorr[2][i], \
            Vspeccorr[3][i], Vspeccorr[4][i], Vspeccorr[5][i])#, Vspeccorr[6][i])   
            Vcorrfile.write(outline)
    else:
        for i in range(4096):
            outline = " %4d %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n" \
               % (i, Ispeccorr[0][i], Ispeccorr[1][i], Ispeccorr[2][i], \
                Ispeccorr[3][i], Ispeccorr[4][i], Ispeccorr[5][i], Ispeccorr[6][i])   
            Icorrfile.write(outline)
    
        for i in range(4096):
            outline = " %4d %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n" \
               % (i, Qspeccorr[0][i], Qspeccorr[1][i], Qspeccorr[2][i], \
                Qspeccorr[3][i], Qspeccorr[4][i], Qspeccorr[5][i], Qspeccorr[6][i])   
            Qcorrfile.write(outline)
        for i in range(4096):
            outline = " %4d %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f \n" \
               % (i, Uspeccorr[0][i], Uspeccorr[1][i], Uspeccorr[2][i], \
                Uspeccorr[3][i], Uspeccorr[4][i], Uspeccorr[5][i], Uspeccorr[6][i])   
            Ucorrfile.write(outline)
        for i in range(4096):
            outline = " %4d %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f \n" \
               % (i, Vspeccorr[0][i], Vspeccorr[1][i], Vspeccorr[2][i], \
                Vspeccorr[3][i], Vspeccorr[4][i], Vspeccorr[5][i], Vspeccorr[6][i])   
            Vcorrfile.write(outline)    
    Icorrfile.close()
    Qcorrfile.close()
    Ucorrfile.close()
    Vcorrfile.close()
    
# plot and fit gain versus declination for each beam
x = array(arange(-5,40,0.1))
alldecs = array(alldecs)
alldecs6 = array(alldecs6)
allKpJy.append(array(allKpJ0))
allKpJy.append(array(allKpJ1))
allKpJy.append(array(allKpJ2))
allKpJy.append(array(allKpJ3))
allKpJy.append(array(allKpJ4))
allKpJy.append(array(allKpJ5))
allKpJy.append(array(allKpJ6))
allrelgain.append(array(relgain1))
allrelgain.append(array(relgain2))
allrelgain.append(array(relgain3))
allrelgain.append(array(relgain4))
allrelgain.append(array(relgain5))
allrelgain.append(array(relgain6))

print "Fit absolute gains"
for i in range(7):
    if(i == 6):
        ax,pconv = curve_fit(func2,alldecs6,allKpJy[i])
        fitx = func2(x,ax[0],ax[1],ax[2])
    else:
        ax,pconv = curve_fit(func4,alldecs,allKpJy[i])
        fitx = func4(x,ax[0],ax[1],ax[2],ax[3],ax[4])
    print ax
    plt.figure(figsize=(10,6))
    plt.title('Beam'+str(i)+' gain. ALFA Angle = 60')
    plt.ylim([5.5,11.5])
    plt.xlabel('Declination (degrees)')
    plt.ylabel('Gain (K/Jy)')
    if(i == 6):
        plt.plot(alldecs6,allKpJy[i],'bo')
    else:
        plt.plot(alldecs,allKpJy[i],'bo')
    plt.plot(x,fitx,'-')
    plt.savefig('Beam'+str(i)+'_gain')
    plt.close()
print "Fit relative gains"    
for i in range(6):
    if (i==5):
        ax,pconv = curve_fit(func2,alldecs6,allrelgain[i])
        fitx = func2(x,ax[0],ax[1],ax[2])
    else:
        ax,pconv = curve_fit(func4,alldecs,allrelgain[i])
        fitx = func4(x,ax[0],ax[1],ax[2],ax[3],ax[4])
    print ax
    plt.figure(figsize=(10,6))
    plt.title('Beam'+str(i+1)+' gain. ALFA Angle = 60')
    plt.ylim([0.5,1.1])
    plt.xlabel('Declination (degrees)')
    plt.ylabel('Gain (K/Jy)')
    if(i == 5):
        plt.plot(alldecs6,allrelgain[i],'bo')
    else:
        plt.plot(alldecs,allrelgain[i],'bo')
    plt.plot(x,fitx,'-')
    plt.savefig('Beam'+str(i+1)+'_relgain')
    plt.close()
   


print "Finished"

