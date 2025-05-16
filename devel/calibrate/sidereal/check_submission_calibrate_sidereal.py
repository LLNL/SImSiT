import yaml
import os

def check_submission(args):
    # Let yaml raise exception here
    print(f"Opening file: {args.infile}")
    with open(args.infile, 'r') as f:
        docs = yaml.safe_load_all(f)
        docs = [doc for doc in docs]

    # First section of yaml file is metadata. Don't include in submission count.
    print(f"Found {len(docs)-1} submissions.")
    # Checking metadata
    assert isinstance(docs[0]['branch'], str)
    assert isinstance(docs[0]['competitor_name'], str)
    assert isinstance(docs[0]['display_true_name'], bool)

    # Do the format check first, since it's fast. Starting after metadata. 
    for doc in docs[1:]:
        # Let python raise error if key not found.
        assert isinstance(doc['file'], str)
        for sat in doc['sats']:
            for var in ['dec0', 'dec1', 'flux', 'mag', 'ra0', 'ra1', 'x0', 'x1', 'y0', 'y1']:
                assert isinstance(sat[var], float)
        if 'stars' in doc:
            for star in doc['stars']:
                for var in ['dec', 'flux', 'mag', 'ra', 'x', 'y']:
                    assert isinstance(star[var], float)

    # Now loop through again if plots requested.
    if args.plot:
        import matplotlib.pyplot as plt
        import astropy.io.fits as fits
        for doc in docs:
            # Again let python raise error if key not found.
            hdul = fits.open(os.path.join(args.public_directory, doc['file']))
            image = hdul[0].data
            ny, nx = image.shape
            kw = {'s':100, 'facecolors':'none'}
            fig, ax = plt.subplots(figsize=(10, 10))
            ax.imshow(image, vmin=0, vmax=200)
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            for sat in doc['sats']:
                ax.scatter([sat['dec0'], sat['dec1']], [sat['ra0'], sat['ra1']], edgecolors='r', **kw)
            if 'stars' in doc:
                for star in doc['stars']:
                    ax.scatter([star['dec']], [star['ra']], edgecolors='c', **kw)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            plt.show()


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("infile")
    parser.add_argument("--plot", action='store_true')
    parser.add_argument("--public_directory", default="public")
    args = parser.parse_args()
    check_submission(args)
