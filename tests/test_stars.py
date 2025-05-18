import numpy as np
import galsim
import satist.tools


def test_mock_star_catalog():
    catalog = satist.catalog.MockStarCatalog()
    rng = np.random.default_rng(57721)
    for _ in range(10):
        z = rng.uniform(-1, 1)
        r = np.sqrt(1-z**2)
        ph = rng.uniform(0, 2*np.pi)
        x = r*np.cos(ph)
        y = r*np.sin(ph)
        center = galsim.CelestialCoord.from_xyz(x, y, z)
        time = None  # Not needed for mock catalog

        stars1 = catalog.get_stars(center, center, 0.1*galsim.degrees, time)
        # Make a second catalog with a displaced center, but still large enough
        # to encompass the original catalog.
        center2 = center.greatCirclePoint(
            galsim.CelestialCoord.from_xyz(0, 0, 1),
            0.09*galsim.degrees
        )
        stars2 = catalog.get_stars(center2, center2, 0.3*galsim.degrees, time)
        assert np.all(np.isin(stars1, stars2))


if __name__ == "__main__":
    test_mock_star_catalog()