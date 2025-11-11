"""
Main entry point for the alleleTools command-line interface.

This module provides the main CLI interface for alleleTools. The CLI is
organized into subcommands for different categories of operations.
"""

import argparse

from . import format, plot


def main():
    parser = argparse.ArgumentParser(
        prog="alleleTools",
        description="Toolset for handling genetic alleles",
    )

    subparsers = parser.add_subparsers(
        dest="command", help="Available commands", required=True
    )

    # Register subcommands
    format.setup_parser(subparsers)
    plot.setup_parser(subparsers)

    args = parser.parse_args()

    if hasattr(args, "func") and args.func is not None:
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
