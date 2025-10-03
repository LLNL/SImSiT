"""
Star catalog interfaces for optical observations of satellites.
"""

import os
from abc import ABC, abstractmethod

import numpy as np

import astropy.units as u
import lsst.sphgeom as sphgeom
from astropy.table import Table, vstack
from ssapy.utils import catalog_to_apparent
from .photometry import get_gaia_magnitude_histogram, vega_to_ab_offset

def _coord_to_UV3d(coord):
    """Convert a coordinate to an LSST spherical geometry UnitVector3d.
    
    Parameters
    ----------
    coord : galsim.CelestialCoord
        Coordinate to convert
        
    Returns
    -------
    sphgeom.UnitVector3d
        Unit vector representation of the coordinate
    """
    return sphgeom.UnitVector3d(*coord.get_xyz())


def _circle_to_poly(circle, nvertex=100):
    """Convert a spherical circle to a convex polygon approximation.
    
    Parameters
    ----------
    circle : sphgeom.Circle
        Spherical circle to convert
    nvertex : int, optional
        Number of vertices to use for polygon approximation (default: 100)
        
    Returns
    -------
    sphgeom.ConvexPolygon
        Polygon approximation of the circle
    """
    center = circle.getCenter()
    # Start with point towards the east and `opening-angle` away.
    point = center.rotatedAround(
        sphgeom.UnitVector3d.northFrom(center),
        circle.getOpeningAngle()
    )
    # Now rotate that point around circle center, mapping out a ~circle.
    th = np.linspace(0, 2*np.pi, nvertex)
    points = [
        point.rotatedAround(center, sphgeom.Angle.fromRadians(th_))
        for th_ in th
    ]
    return sphgeom.ConvexPolygon(points)


class StarCatalog(ABC):
    """Abstract base class for star catalogs using HTM pixelization.
    
    Parameters
    ----------
    level : int, optional
        HTM pixelization level (default: 7)
    """
    def __init__(self, level=7):
        self.htm = sphgeom.HtmPixelization(level)

    def get_stars(self, coord0, coord1, radius, time):
        """Get stars within a field of view, accounting for motion during exposure.
        
        Parameters
        ----------
        coord0, coord1 : galsim.CelestialCoord
            Field of view center at beginning/end of exposure.
        radius : galsim.Angle
            Field of view radius
        time : astropy.time.Time
            Time of observation (to account for proper motion, parallax,
            aberration)

        Returns
        -------
        stars : astropy.table.Table
            Table with columns
                - ra (deg)
                - dec (deg)
                - i_mag
        """
        c0 = sphgeom.Circle(
            _coord_to_UV3d(coord0),
            sphgeom.Angle.fromDegrees(radius.deg)
        )
        c1 = sphgeom.Circle(
            _coord_to_UV3d(coord1),
            sphgeom.Angle.fromDegrees(radius.deg)
        )
        poly0 = _circle_to_poly(c0)
        poly1 = _circle_to_poly(c1)
        poly = sphgeom.ConvexPolygon(poly0.getVertices() + poly1.getVertices())

        htm = sphgeom.HtmPixelization(7)
        ranges = htm.envelope(poly)
        table_list = []
        for start, end in ranges:
            for idx in range(start, end):
                table_list.append(self._get_trixel_stars(idx, time))
        out = vstack(table_list)
        w = poly.contains(out['ra_rad'], out['dec_rad'])
        out = out[w]
        return out

    @abstractmethod
    def _get_trixel_stars(self, idx, time):
        """Get stars from a single HTM trixel.
        
        Parameters
        ----------
        idx : int
            HTM trixel index
        time : astropy.time.Time
            Time of observation
            
        Returns
        -------
        astropy.table.Table
            Table of stars in the trixel
        """
        pass


class MockStarCatalog(StarCatalog):
    """A mock star catalog created on-the-fly in memory, but consistently (same
    stars in the same positions regardless of query parameters).

    Parameters
    ----------
    level : int, optional
        HTM level for pixelization (default: 7)
    seed : int, optional
        Random seed for initialization (default: 123)

    Notes
    -----
    Catalog is created on the fly in each HTM pixel (trixel), based on that
    trixel's id and the given random seed.  So for a given level and seed, the
    star catalog is fixed.  With a different level or seed, the catalog will
    change.
    """
    def __init__(self, level=7, seed=123):
        super().__init__(level=level)
        self.seed = seed
        # i magnitude histogram from limited area GAIA trixels
        self.imag_hist = np.array([
            129,   168,   180,   170,   191,   213,   233,   261,   249,
            300,   338,   350,   393,   413,   437,   476,   502,   504,
            596,   633,   687,   717,   753,   850,   898,   980,  1027,
            1122,  1110,  1230,  1347,  1382,  1418,  1580,  1687,  1774,
            1818,  1972,  2058,  2173,  2275,  2419,  2544,  2679,  2905,
            3062,  3135,  3286,  3469,  3619,  3748,  3894,  4079,  4334,
            4599,  4736,  4898,  5114,  5577,  5695,  5772,  5992,  6346,
            6562,  6716,  7072,  7293,  7516,  7833,  8116,  8161,  8800,
            9135,  9282,  9633, 10074, 10294, 10719, 11378, 11609, 12084,
            12509, 12647, 13454, 14090, 14407, 15194, 15699, 16437, 16951,
            17812, 18275, 19097, 20012, 20865, 21666, 22424, 23508, 24474,
            25621, 26306, 27156, 28272, 28551, 29535, 30510, 31067, 32346,
            33013, 33201
        ])
        self.bins = np.linspace(10, 21, 111)
        self.density = np.sum(self.imag_hist) / 0.102219849198379 # N / sr

        self.magnitude_hists = get_gaia_magnitude_histogram()
       

    def _get_trixel_stars(self, idx, time):
        """Generate mock stars for a given HTM trixel.
        
        Generates a deterministic random star catalog based on the trixel index
        and seed. Stars are uniformly distributed within the trixel with magnitudes
        sampled from GAIA magnitude histograms.
        
        Parameters
        ----------
        idx : int
            HTM trixel index
        time : astropy.time.Time
            Time of observation (unused in mock catalog)
            
        Returns
        -------
        astropy.table.Table
            Table with columns:
                - ra_rad : RA in radians
                - dec_rad : Dec in radians
                - ra : RA in degrees
                - dec : Dec in degrees
                - i_mag : i-band magnitude
                - {filter}_mag : Magnitude in GAIA filters
        """
        # Generate random deterministic star catalog anywhere on the sky on the
        # fly.
        trixel = self.htm.triangle(idx)
        circle = trixel.getBoundingCircle()
        area = circle.getArea()

        rng = np.random.default_rng(idx+self.seed)
        N = int(rng.poisson(self.density*area))

        # populate spherical cap
        h = area/(2*np.pi)
        z = rng.uniform(1-h, 1, size=N)
        ph = rng.uniform(0, 2*np.pi, size=N)
        r = np.sqrt(1-z**2)
        x = r * np.cos(ph)
        y = r * np.sin(ph)
        r = np.array([x, y, z])

        # Now rotate to desired center
        k = np.cross([0,0,1], circle.getCenter())
        k /= np.sqrt(np.sum(k**2))
        th = np.pi/2 - np.arcsin(circle.getCenter()[2])
        sth, cth = np.sin(th), np.cos(th)
        # Rodriguez rotation formula
        r = cth*r + sth*(np.cross(k, r.T).T) + (1-cth)*np.dot(k, r)*k[:,None]
        x, y, z = r

        ra = np.arctan2(y, x)
        dec = np.arcsin(z)
        w = trixel.contains(ra, dec)

        table = Table()
        table['ra_rad'] = ra[w]
        table['dec_rad'] = dec[w]

        # Roughly sample the histogram in __init__.
        # Uniform probability in each fixed-sized bin
        indices = rng.choice(
            len(self.imag_hist),
            p=self.imag_hist/np.sum(self.imag_hist),
            size=N
        )
        i_mag = self.bins[indices] + rng.uniform(0, 0.1, size=N)
        table['i_mag'] = i_mag[w]

        for gaia_filter in self.magnitude_hists:
        
            mag_hist = self.magnitude_hists[gaia_filter]['values']
            bins = self.magnitude_hists[gaia_filter]['bins']
            indices = rng.choice(
                len(mag_hist),
                p=mag_hist/np.sum(mag_hist),
                size=N
            )
            mag = bins[indices] + rng.uniform(0, 0.1, size=N)
            table[f'{filter}_mag'] = mag[w] + vega_to_ab_offset(gaia_filter)


        table['ra'] = np.rad2deg(table['ra_rad'])
        table['dec'] = np.rad2deg(table['dec_rad'])

        return table


class GaiaStarCatalog(StarCatalog):
    """Get star catalog from GAIA and convert from ICRF to apparent coordinates.
    
    Reads GAIA catalog data from HTM-pixelized FITS files and transforms
    coordinates to apparent positions accounting for proper motion, parallax,
    and aberration.

    Parameters
    ----------
    gaia_dir : str
        Directory containing GAIA catalog trixel FITS files (named {idx}.fits)
    level : int, optional
        HTM pixelization level (default: 7)
    """
    def __init__(self, gaia_dir, level=7):
        super().__init__(level=level)
        self.gaia_dir = gaia_dir

    def _get_trixel_stars(self, idx, time):
        """Load and transform GAIA stars from a single HTM trixel.
        
        Reads GAIA catalog data from disk and applies proper motion, parallax,
        and aberration corrections to transform from ICRF to apparent coordinates.
        
        Parameters
        ----------
        idx : int
            HTM trixel index
        time : astropy.time.Time
            Time of observation for coordinate transformations
            
        Returns
        -------
        astropy.table.Table
            Table with columns:
                - ra_rad : Apparent RA in radians
                - dec_rad : Apparent Dec in radians
                - ra : Apparent RA in degrees
                - dec : Apparent Dec in degrees
                - i_mag : i-band magnitude (currently G-band flux)
        """
        file = os.path.join(self.gaia_dir, f"{idx}.fits")
        gaia_data = Table.read(file)
        table = Table()
        # Apply transformations for proper motion, parallax, and aberration.
        # The pm and px transforms are tiny, but easy to include.  For
        # aberration, we are transforming into the center-of-earth velocity
        # frame, which is the frame we use to calculate orbits.  We could
        # also include diurnal aberration and refraction here, but since these
        # equally affect stars and satellites, we simply leave them out and
        # don't calculate them for satellites either.
        table['ra_rad'], table['dec_rad'] = catalog_to_apparent(
            gaia_data['coord_ra'].to(u.rad).value,
            gaia_data['coord_dec'].to(u.rad).value,
            time,
            observer=None,  # we don't want to include diurnal aberration
            pmra=gaia_data['pm_ra'].to(u.mas).value,
            pmdec=gaia_data['pm_dec'].to(u.mas).value,
            parallax=gaia_data['parallax'].to(u.arcsec).value,
        )
        gaia_data = gaia_data[gaia_data['phot_g_mean_flux'] > 0.0]
        # TODO: change line below for multiple filters
        # g_coefficients = []
        # r_coefficients = []
        # i_coefficients = []
        # TODO: logic for which filter
        # if filter == 'i':
        #    coefficients = i_coefficients
        #    filter_name = 'i_mag'
        # deltaG = 
        # table[filter_name] = -(coefficients[0] - (coefficients[1] * deltaG))
        table['i_mag'] = gaia_data['phot_g_mean_flux'].to(u.ABmag).value
        table['ra'] = np.rad2deg(table['ra_rad'])
        table['dec'] = np.rad2deg(table['dec_rad'])
        return table