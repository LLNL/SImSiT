"""
Telescope and detector instrument models for satellite observation simulations.
"""

import galsim
import numpy as np

from .wcs import radialWCS


class Instrument:
    """Instrument model for telescope and detector properties.

    Parameters
    ----------
    image_shape : 2-tuple of int
        Image dimensions (nx, ny) in pixels
    gain : float
        Detector gain in electrons per ADU
    read_noise : float
        Read noise in electrons
    pixel_scale : float
        Physical pixel size in microns
    aperture : float
        Telescope aperture diameter in meters
    obscuration : float, optional
        Linear fractional obscuration (default: 0.0)
    distortion : dict, optional
        Optical distortion model with fields:
            'th' : array
                Field angle in degrees
            'dthdr' : array
                Plate scale in arcsec/micron
    vignetting : dict, optional
        Vignetting model with fields:
            'th' : array
                Field angle in degrees
            'unvig' : array
                Surviving fraction of photons
    """
    def __init__(
        self,
        *,
        image_shape,
        gain,
        read_noise,
        pixel_scale,
        aperture,
        obscuration=0.0,
        distortion=None,
        vignetting=None
    ):
        self.image_shape = image_shape
        self.gain = gain
        self.read_noise = read_noise
        self.pixel_scale = pixel_scale
        self.aperture = aperture
        self.obscuration = obscuration
        self.distortion = distortion
        self.vignetting = vignetting

        if vignetting is not None:
            from scipy.interpolate import interp1d
            self.vigfun = interp1d(
                vignetting['th'],
                vignetting['unvig'],
                kind='cubic'
            )

    @staticmethod
    def fromConfig(config):
        """Create an Instrument instance from a configuration dictionary.
        
        Parameters
        ----------
        config : dict
            Configuration dictionary with required keys:
                - image_shape : 2-tuple of int
                - gain : float
                - read_noise : float
                - pixel_scale : float
                - aperture : float
            Optional keys:
                - obscuration : float (default: 0.0)
                - distortion : dict
                - vignetting : dict
                
        Returns
        -------
        Instrument
            Configured instrument instance
        """
        return Instrument(
            image_shape = config['image_shape'],
            gain = config['gain'],
            read_noise = config['read_noise'],
            pixel_scale = config['pixel_scale'],
            aperture = config['aperture'],
            obscuration = config.get('obscuration', 0.0),
            distortion = config.get('distortion', None),
            vignetting = config.get('vignetting', None)
        )

    def init_image(self, *, sky_phot, exptime):
        """Initialize an image with sky background.

        Parameters
        ----------
        sky_phot : float
            Sky level in photons / arcsec^2 / sec
        exptime : float
            Exposure time in seconds

        Returns
        -------
        image : galsim.Image
            Image with sky background added
        """
        nx, ny = self.image_shape
        bounds = galsim.BoundsI(-nx//2, nx//2-1, -ny//2, ny//2-1)
        image = galsim.Image(bounds)
        # Get a mock wcs by asserting arbitrary boresight, rotation.  Only
        # actually care about relative pixel sizes here.
        boresight = galsim.CelestialCoord(0*galsim.degrees, 0*galsim.degrees)
        rot_sky_pos = 0*galsim.degrees
        mock_wcs = self.get_wcs(boresight, rot_sky_pos)
        # Sky background
        mock_wcs.makeSkyImage(image, sky_phot*exptime)
        return image

    @galsim.utilities.lazy_property
    def field_radius(self):
        """Calculate the field radius from boresight to corner.
        
        Returns
        -------
        galsim.Angle
            Angular distance from boresight to image corner
        """
        boresight = galsim.CelestialCoord(0*galsim.degrees, 0*galsim.degrees)
        rot_sky_pos = 0*galsim.degrees
        mock_wcs = self.get_wcs(boresight, rot_sky_pos)
        nx, ny = self.image_shape
        coord = mock_wcs.toWorld(galsim.PositionD(nx/2, ny/2))
        return coord.distanceTo(boresight)

    def apply_vignetting(self, image):
        """Apply vignetting correction to an image.
        
        Modifies the image in-place by multiplying each pixel by the vignetting
        function based on radial distance from the center.
        
        Parameters
        ----------
        image : galsim.Image
            Image to apply vignetting to (modified in-place)
        """
        # Vignette the background
        # Should really use wcs here to use angular distances, but for now we'll
        # cheat and just use pixel distances.
        boresight = galsim.CelestialCoord(0*galsim.degrees, 0*galsim.degrees)
        rot_sky_pos = 0*galsim.degrees
        mock_wcs = self.get_wcs(boresight, rot_sky_pos)
        xx, yy = image.get_pixel_centers()
        rr = np.hypot(xx, yy)  # dist from center in pixels
        rr *= np.sqrt(mock_wcs.pixelArea(galsim.PositionD(0, 0)))  # -> arcsec
        image.array[:] *= self.vigfun(rr/3600)

    def apply_noise(self, image, sky_phot, rng):
        """Add CCD noise (Poisson, read noise) to an image.
        
        Parameters
        ----------
        image : galsim.Image
            Image to add noise to (modified in-place)
        sky_phot : float
            Sky level in photons / arcsec^2 / sec
        rng : numpy.random.Generator
            Random number generator for noise generation
        """
        gsrng = galsim.BaseDeviate(rng.bit_generator.random_raw() % (2**63))
        noise = galsim.CCDNoise(
            gsrng,
            sky_level=sky_phot*self.pix_size**2,
            gain=self.gain,
            read_noise=self.read_noise
        )
        image.addNoise(noise)

    def get_wcs(self, boresight, rot_sky_pos):
        """Create a WCS object for the instrument.
        
        Parameters
        ----------
        boresight : galsim.CelestialCoord
            Telescope pointing direction (field center)
        rot_sky_pos : galsim.Angle
            Position angle of detector "up" direction with respect to North

        Returns
        -------
        wcs : galsim.GSFitsWCS
            World coordinate system transformation
        """
        return radialWCS(
            self.distortion['th'],
            np.array(self.distortion['dthdr'])*self.pixel_scale,
            boresight,
            rot_sky_pos,
            n=8, order=3
        )

    @galsim.utilities.lazy_property
    def pix_size(self):
        """Calculate the mean pixel size in arcseconds.
        
        Returns
        -------
        float
            Mean pixel size in arcseconds
        """
        return np.mean(self.distortion['dthdr']) * self.pixel_scale

    def compute_LSST_scaled_zp(self):
        """Compute the photometric zero point scaled from LSST values.
        
        Scales the LSST i-band zero point based on the ratio of collecting areas
        between this instrument and LSST, accounting for central obscuration.
        
        Returns
        -------
        float
            Photometric zero point in AB magnitudes, where an object of magnitude
            zp produces 1 photon/sec
            
        Notes
        -----
        Based on LSST i-band values from:
        https://github.com/LSSTDESC/WeakLensingDeblending/blob/master/descwl/survey.py
        This scaling only works for LSST filters; SWIR bands require different treatment.
        """
        # Scale from the LSST i-band values listed here:
        # https://github.com/LSSTDESC/WeakLensingDeblending/blob/master/descwl/survey.py
        LSST_area = 32.4  # meters^2
        # TODO: edit line below for other bands. This scaling only works for LSST filters, going into SWIR we won't have this.
        LSST_ZP = 32.36  # electrons / sec at i_AB = 24
        # Smaller telescope will receive proportionally fewer photons for the same exptime
        radius = self.aperture*0.5
        area = np.pi*(radius**2*(1-self.obscuration**2))
        ZP = LSST_ZP * area/LSST_area

        # Convert `ZP` to more traditional `zp` where object of mag zp produces 1 phot / sec
        # Solve
        #     ZP 10^(-0.4 (i - 24)) = 10^(-0.4 (i - zp))
        #     log10(ZP) - 0.4 i + 0.4 24 = -0.4 i + 0.4 zp
        #     zp = log10(ZP)/0.4 + 24
        return  np.log10(ZP)/0.4 + 24

    def streak_snr(self, *, nphot, length, psf_fwhm, sky_phot):
        """Calculate signal-to-noise ratio for a streak.
        
        Computes SNR for a linear streak using an effective aperture that accounts
        for the PSF size and streak length.
        
        Parameters
        ----------
        nphot : array_like
            Number of photons in the streak
        length : float
            Streak length in arcseconds
        psf_fwhm : float
            PSF full-width at half-maximum in arcseconds
        sky_phot : float
            Sky level in photons / arcsec^2 / sec

        Returns
        -------
        snr : array_like
            Signal-to-noise ratio for each streak
        """
        neff = 2.266 * (psf_fwhm / self.pix_size)**2
        reff = np.sqrt(neff/np.pi)
        # effective area is rectangle + 2 semi-circular endcaps
        aeff = neff+reff*length
        var = nphot/self.gain
        var += (sky_phot*self.pix_size**2/self.gain + self.read_noise**2)*aeff
        return nphot / np.sqrt(var)