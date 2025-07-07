from . import graph_pathogens, graph_phewas


def setup_parser(subparsers):
    parser = subparsers.add_parser(
        "interpret",
        help="Get useful information about HLA alleles.",
        description="This commands provide some useful tools"
                    "for data interpretation. e.i. pleiotropy reports",
    )
    interpret_parser = parser.add_subparsers(
        dest="interpret_type",
        help="Available commands:",
        required=True,
    )

    graph_phewas.setup_parser(interpret_parser)
    graph_pathogens.setup_parser(interpret_parser)

    return parser
