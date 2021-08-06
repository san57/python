import argparse

from pycif.utils.classes.setup import Setup


def main():
    """Main entry point to pyCIF.
    The only required argument is the defintion Yaml file

    """

    # Parses the arguments given to the script
    # Use the following command to see the script helper:
    # $python -m pycif -h
    cmd_parser = argparse.ArgumentParser(
        description="Runs pyCIF from a given Yaml file"
    )
    cmd_parser.add_argument(
        "def_file", type=str, help="path to the Yaml configuration file"
    )
    cmd_parser.add_argument("--debug", help="debug mode", action="store_true")

    args = vars(cmd_parser.parse_args())

    # Run pyCIF
    Setup.run_simu(args)


if __name__ == "__main__":
    main()
