import yaml
import sys 
#sys.path.insert(1, "./")

from satist import simulate


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument(
        "config",
        type=str
    )
    args = parser.parse_args()

    config = yaml.safe_load(open(args.config, 'r'))
    simulate.simulate(config)
    simulate.make_sky_flat(args.config)
