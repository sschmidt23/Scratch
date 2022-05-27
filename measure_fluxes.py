import numpy as np
import scipy.interpolate as sciint

class Spectrum(object):
    """
    Class to take a spectrum and compute the flux in a set of filters
    Inputs:
    -------
    wl: np like float array
       Array of wavelength values, assume in Angstroms
    flux: np like float array
       Array of fluxes, assume in f_lambda units (as opposed to f_nu)
    filterlist: list of strings
       An array of the filenames for the filters to be read in
    """
    def __init__(self,wl,flux,filterlist):
        self.wl = wl
        self.flux = flux
        self.flux_fnu = flux*wl*wl
        self.filterlist = filterlist
        filter_wls,filter_trans = self.ReadFilters()
        self.filter_wls = filter_wls
        self.filter_trans = filter_trans
        
    def ReadFilters(self):
        """
        Reads in filters and returns a dictionary of wls and a dictionary of
        transmissions for the filters
        Returns:
        --------
        filter_wl: dictionary
           dictionary of float arrays of wavelengths for each transmission curve
        filter_trans: dictionary
           dictionary of float arrays of transmissions for each transmission 
           curve
        """
        filterlist = self.filterlist
        num_filts = len(filterlist)
        filter_wl = {}
        filter_trans = {}
        for filt in filterlist:
            filtname = filt[:-4]
            print (filtname)
            data = np.loadtxt(filt)
            indiv_wl = data[:,0]
            indiv_trans = data[:,1]
            #norm = np.sum(indiv_wl*indiv_trans)
            #indiv_trans /= norm
            filter_wl[filtname] = indiv_wl
            filter_trans[filtname] = indiv_trans
        return filter_wl,filter_trans

    
    def PrintFilters(self):
        for filt in self.filterlist:
            print (filt)
        return
    
    def RedshiftSpectrum(self,z,interpolate_flag = True):
        """
        function to redshift spectrum to redshift z.  If interpolate_flag is
        True, re-sample the spectrum on to the original wavelength grid
        input:
        ------
        z: float
           redshift to shift spectrum to
        interpolate_flag: bool
           interpolate fluxes back to original wl grid if True
        Returns:
        --------
        interp_wl: array
           the wavelengths of the redshifted spectrum
        flux:
           the same fluxes in self.flux, (just for consistency)
        """
        self.z = z
        new_wl = self.wl*(1.+z)
        if interpolate_flag:
            interp_func = sciint.interp1d(new_wl,self.flux,kind = 'linear',
                                        bounds_error=False,fill_value=0.0)

            interp_fl = [interp_func(x) for x in self.wl]
            self.z_wl = self.wl
            self.z_flux = interp_fl
            return self.wl,interp_fl
        else:
            self.z_wl = new_wl
            self.z_flux = self.flux
            return new_wl,self.flux

    def CalculateRedshiftedFilterFlux(self,xfilter):
        """
        calculate the flux * transmission in given filter
        inputs:
        -------
        xfilter: specific filter, which must be in filterlist
        Returns:
        --------
        filt_flux: float
        flux*transmission for the specific filter
        """
        filtwl = self.filter_wls[xfilter]
        filttrans = self.filter_trans[xfilter]
        #need to get spectrum on same grid as filter
        #calculate f_nu for redshifted spectrum
        z_fl_fnu = self.z_flux*self.z_wl*self.z_wl
        interp_func = sciint.interp1d(self.z_wl,z_fl_fnu,kind = 'linear',
                                      bounds_error=False,fill_value=0.0)
        norm = np.trapz(filttrans,filtwl)
        interp_flux = [interp_func(x) for x in filtwl]
        filt_flux = np.trapz(filttrans*interp_flux,filtwl)/norm
        #lp = np.sqrt(norm/np.trapz(filttrans/(filtwl*filtwl),filtwl))
        #filt_flux = filt_flux_lam*lp**2
        return filt_flux
