import numpy as np
import measure_fluxes as meas_flux
import matplotlib.pyplot as plt
import sys,os

def main():

    basepath = "./"
    filebase = "_dc2_ugsorted_v4.sed"
    filterpath = "./"
    filterlist = ["DC2LSST_u","DC2LSST_g","DC2LSST_r","DC2LSST_i","DC2LSST_z","DC2LSST_y"]

    filterlistx = ["DC2LSST_u.res","DC2LSST_g.res","DC2LSST_r.res","DC2LSST_i.res","DC2LSST_z.res","DC2LSST_y.res"]

    z_arr = np.arange(0.0,3.051,0.1)
    #flux_array = np.zeros([len(filterlist),len(z_arr)])
    color_array = np.zeros([len(filterlist)-1,len(z_arr)])

    outfp = open("cosmodc2template_wbc03_v4_colors.out","w")
    outfp.write("#templatenum redshift colors\n")
    for i in range(229):
        tmpname = "%d"%i+filebase
        fullpath = os.path.join(basepath,tmpname)
        data = np.loadtxt(fullpath,skiprows=1)
        wl = data[:,0]
        fl = data[:,1]

        spectrum = meas_flux.Spectrum(wl,fl,filterlistx)
        print ("loaded spectrum %d"%i)
        #print ("following filters are loaded:")
        #spectrum.PrintFilters()

        for j,zz in enumerate(z_arr):
            print zz
            outfp.write("%d %.2f "%(i, zz))
            flux_array = np.zeros([len(filterlist)])
            spec_wl,spec_fl = spectrum.RedshiftSpectrum(zz,
                                                        interpolate_flag=True)
            for k,filt in enumerate(filterlist):
                flux_array[k] = spectrum.CalculateRedshiftedFilterFlux(filt)
                
            for k in range(len(filterlist)-1):
                color_array[k,j] = -2.5*np.log10(flux_array[k]/flux_array[k+1])
                outfp.write("%.4f "%(color_array[k,j]))
            outfp.write("\n")
    outfp.close()
  

  
#    fig = plt.figure(figsize=(8,12))
#    ax = plt.subplot(411)
#    plt.plot(color_array[0,:],color_array[1,:],marker='.',s=10,c='r')
#    plt.xlabel("u-g",fontsize=12)
#    plt.ylabel("g-r",fontsize=12)
#    plt.show()

    
if __name__=="__main__":
    main()
