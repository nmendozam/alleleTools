import argparse
from . import convert


def main():
    parser = argparse.ArgumentParser(
            prog='alleleTools',
            description='Toolset for handling genetic alleles',
            )

    subparsers = parser.add_subparsers(
            dest='command',
            help='Available commands',
            required=True
            )

    convert.setup_parser(subparsers)

    args = parser.parse_args()

    if hasattr(args, 'func') and args.func is not None:
        args.func(args)
    else:
        parser.print_help()

if __name__ == '__main__':
    main()
