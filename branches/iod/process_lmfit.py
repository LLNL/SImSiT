#
# LLNL method for IOD
#
import numpy as np
import astropy.units as u
from astropy.table import QTable, Table
from astropy.time import Time
import ssa
import functools
import pdb


def iod_lmfit(track, propagator=ssa.SGP4Propagator()):
    """
    An IOD method based on LM optimization in the LLNL ssa code

    Parameters
    ----------
    track : astropy.Table
	propagator : ssa.Propagator derived class

    Returns
    -------
    rv : array_like (6,)
        r and v at epoch (meters and meters / second)
    epoch : float
        Epoch time in gps seconds
    """
    # print(observer.name)
    # rStation, _ = observer.getRV(track['time'])

    arc = QTable()
    arc['ra'] = track['ra'] * u.deg
    arc['dec'] = track['dec'] * u.deg
    arc['time'] = Time(track['time'], format='gps')
    arc['rStation_GCRF'] = np.column_stack((track['px'], track['py'], track['pz'])) * u.m
    arc['sigma'] = 10.0 * u.arcsec

    # this ignores precession, etc., but I want to get something in here.
    # also assumes that observers are on the earth
    # note that velocity is only important for things like aberration.
    l, b = ssa.utils.unit2lb(arc['rStation_GCRF'].value)
    up = np.array([0, 0, 1])
    rotationrate = ssa.constants.WGS84_EARTH_OMEGA/u.s
    vv = np.cross(up[None, :]*rotationrate, arc['rStation_GCRF'], axis=1)
    arc['vStation_GCRF'] = vv

    priors = [ssa.rvsampler.EquinoctialExEyPrior(0.3)]
    # circular guess needs the first two things in arc to be from the same
    # site relatively close in time to one another.
    # this is not always the case in this data set.
    # ugly hacks follow.
    rsite = np.sqrt(np.sum(arc['rStation_GCRF'].to(u.km).value**2, axis=1)).round(3)
    rsiteu = np.unique(rsite)
    nmax = 0
    for rr in rsiteu:
        nn = np.sum(rr == rsite)
        if nn > nmax:
            nmax = nn
            rsite_mostcommon = rr
    m = rsite == rsite_mostcommon
    dt = np.diff(arc['time'].gps[m])
    bestind = np.argmin(dt + (dt <= 0)*1e9)

    initRV, epoch = ssa.circular_guess(arc[m][bestind:bestind+2])
    probfn = ssa.RVProbability(arc, epoch, propagator=propagator,
                               priors=priors)
    initEl = ssa.Orbit(r=initRV[:3], v=initRV[3:6], t=epoch).equinoctialElements
    lmo = ssa.rvsampler.SGP4LMOptimizer(probfn, initEl)
    el = lmo.optimize()
    orb = lmo._getOrbit(el)
    epochfinal = np.max(arc['time'].gps)
    probfnfinal = ssa.RVProbability(arc, epochfinal, propagator=propagator,
                                    priors=priors)
    rfinal, vfinal = ssa.rv(orb, epochfinal, propagator=propagator)
    rvfinal = np.concatenate([rfinal, vfinal])
    lmofinal = ssa.rvsampler.LMOptimizer(probfnfinal, rvfinal)
    lmofinal.optimize()
    covar = lmofinal.result.covar
    resid = lmofinal.result.residual
    ok = np.sum(resid**2)/len(resid) < 10

    up = np.array([0, 0, 1])
    rotationrate = ssa.constants.WGS84_EARTH_OMEGA
    prop_times = [5, 15, 30, 60]
    def fictionalrd(param, dt):
        orb = ssa.Orbit(r=param[:, :3], v=param[:, 3:6], t=epochfinal)
        r, v = ssa.compute.rv(orb, epochfinal+dt, propagator=propagator)
        rGround = ssa.utils.normed(r)*ssa.constants.WGS84_EARTH_RADIUS
        vGround = np.cross(up[None, :]*rotationrate, rGround, axis=1)
        ra, dec, _, pmra, pmdec, _ = ssa.compute.rvObsToRaDecRate(
            r, v, rGround, vGround)
        return np.concatenate([[ra], [dec]]).T

    predloc = {}
    for dt in prop_times:
        func = functools.partial(fictionalrd, dt=dt*60)
        fsigma = ssa.utils.sigma_points(func, rvfinal, covar)
        fsigma = np.degrees(fsigma)
        rd = fsigma[0]
        dmean = fsigma[1:, ...] - rd[:, ...]
        cd = np.cos(np.radians(rd[1]))
        dmean[:, 0] = (((dmean[:, 0] + 180) % 360) - 180)*cd
        rdcovar = np.cov(dmean.T, ddof=0)
        predloc[dt] = (rd, rdcovar)

    if not ok:
        return None, None, None, None, None, None
    else:
        return rvfinal, epochfinal, covar, lmofinal, predloc


def main():
    from argparse import ArgumentParser
    from tqdm import tqdm
    import glob
    import os
    import yaml

    parser = ArgumentParser()
    parser.add_argument("--public_directory", type=str, default="public/iod_datasets/")
    parser.add_argument("--outfile", type=str, default="lmfit_iod.yaml")
    args = parser.parse_args()

    files = sorted(glob.glob(os.path.join(args.public_directory,
                                          "SAT_????.csv")))

    docs_with_optional = [{'branch' : 'IOD', 'competitor_name' : 'Competitor A', 'display_true_name' : True}]
    docs = [{'branch' : 'IOD', 'competitor_name' : 'Competitor A', 'display_true_name' : True}]

    propagator = ssa.propagator.SGP4Propagator()
    # propagator = ssa.propagator.KeplerianPropagator()

    for idx, file in enumerate(tqdm(files)):
        track = Table.read(file)
        rv, epoch_gps, covar, _, predloc = iod_lmfit(
            track, propagator=propagator)
        if rv is None:
            continue

        covar = [[float(x) for x in r] for r in covar]
        pred_loc_submit = dict()
        pred_loc_covar_submit = dict()
        for dt in predloc:
            loc = dict()
            loc_covar = dict()
            
            loc['ra'] = float(predloc[dt][0][0])
            loc['dec'] = float(predloc[dt][0][1])
            loc_covar = [[float(x) for x in r] for r in predloc[dt][1]]
            pred_loc_covar_submit[dt] = loc_covar
            pred_loc_submit[dt] = loc

        iod = dict()
        loc = dict()
        iod['rx'] = float(rv[0])
        iod['ry'] = float(rv[1])
        iod['rz'] = float(rv[2])
        iod['t'] = float(epoch_gps)
        iod['vx'] = float(rv[3])
        iod['vy'] = float(rv[4])
        iod['vz'] = float(rv[5])
        loc = pred_loc_submit
        state_vector_covar = covar
        predicted_location_covar = pred_loc_covar_submit
        docs_with_optional.append({'file': os.path.basename(file), 'IOD':[iod], 'predicted_downrange_location':[loc], 'state_vector_covariance':[covar], 'predicted_location_covariance':[pred_loc_covar_submit]})
        docs.append({'file': os.path.basename(file), 'IOD':[iod]})

        if idx == 4:
            with open("public/sample_submission_5_iod.yaml", "w") as f:
                yaml.safe_dump_all(docs_with_optional, f)
            with open("public/sample_submission_5_iod_required.yaml", "w") as f:
                yaml.safe_dump_all(docs, f)

    with open("private/sample_submission.yaml", "w") as f:
        yaml.safe_dump_all(docs_with_optional, f)
    with open("private/sample_submission_required.yaml", "w") as f:
        yaml.safe_dump_all(docs, f)

    with open(args.outfile, "w") as f:
    	yaml.safe_dump_all(docs_with_optional, f)


if __name__ == '__main__':
    main()

