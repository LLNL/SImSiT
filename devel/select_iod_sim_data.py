#!/usr/bin/env python
#
# Parse IOD simulation data to make individual track subsets for validation tests
#
import os.path
import numpy as np
from astropy.table import Table

k_rootdir = './'


def save_tracks(datadir, num_sats=2):
    for satnum in range(num_sats):
        obs = Table.read(os.path.join(datadir, f'{satnum:03}_observation.fits'), format='fits')
        eph = Table.read(os.path.join(datadir, f'{satnum:03}_emphemeris.fits'), format='fits')

        out_head = os.path.join(k_rootdir, 'data/llnl_sim_tracks', f'sat{satnum:03}')
        print("Saving to ", out_head)
        os.makedirs(out_head, exist_ok=True)

        obs_sites = obs.group_by('obs_site')
        for site in obs_sites.groups:
            obs_track = site.group_by('obs_track_id')
            for g in obs_track.groups:
                track = g['obs_site','obs_t','obs_ra','obs_dec']

                track_file_label = f"{g['obs_site'][0]}_{g['obs_track_id'][0]:03d}"
                ### Save ephemeris at future times for truth comparisons
                ###   We check the state at 30-34 minutes downrange of the first obs
                ndx = np.logical_and(eph['eph_t'] > g['obs_t'][0] + 30*60,
                                     eph['eph_t'] < g['obs_t'][0] + 34*60)

                truthdir = os.path.join(k_rootdir, 'truths/llnl_sim_tracks', f'sat{satnum:03}')
                os.makedirs(truthdir, exist_ok=True)
                truthfile = os.path.join(truthdir, track_file_label + ".fits")
                truth = eph[ndx]
                truth.write(truthfile, format='fits')

                timesdir = os.path.join(out_head, 'pred_times')
                os.makedirs(timesdir, exist_ok=True)
                timesfile = os.path.join(timesdir, track_file_label + '.txt')
                pred_times = truth['eph_t']
                np.savetxt(timesfile, pred_times)

                ### Save track subsets of varying lengths for testing
                for num_obs_per_track in range(2, len(track)):
                    outdir = os.path.join(out_head, f'tracks{num_obs_per_track:03}')
                    os.makedirs(outdir, exist_ok=True)
                    outfile = os.path.join(outdir, track_file_label + ".txt")
                    track[0:num_obs_per_track].write(outfile, format='ascii')
    return None

def main():
    datadir = 'data/sim_for_iod'

    ### Save the tracks
    save_tracks(datadir, num_sats=2)


if __name__ == '__main__':
    main()