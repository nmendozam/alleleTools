from . import consensus


def setup_parser(subparsers):
    parser = subparsers.add_parser(
        'genotype',
        help='Useful commands to use when genotyping polymorphic genes.',
        description=
                    "This command has some useful algorithms to use when genotyping polymorphic genes."
    )
    genotype_parser = parser.add_subparsers(
            dest='genotype_type',
            help='Specify the operation to perform.',
            required=True,
    )

    consensus.setup_parser(genotype_parser)

    return parser
