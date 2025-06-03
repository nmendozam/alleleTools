# from . import allele_resolution


def setup_parser(subparsers):
    parser = subparsers.add_parser(
        'refactor',
        help='Refactor allele table data.',
        description=
                    "This command facilitates usual refactoring of allele table data."
    )
    refactor_parser = parser.add_subparsers(
            dest='refactor_type',
            help='Specify the type of refactoring to perform.',
            required=True,
    )

    # allele_resolution.setup_parser(refactor_parser)

    return parser
