from . import vcf2allele


def setup_parser(subparsers):
    parser = subparsers.add_parser(
        'convert',
        help='Convert allele data between different formats.',
        description='Detailed help for the convert command:\n\n'
                    'This command facilitates the conversion of allele data from one\n'
                    'structured format to another (e.g., VCF to CSV, custom format to standard).\n'
                    'It supports various input/output formats and allows for optional\n'
                    'validation during conversion.'
    )
    convert_parser = parser.add_subparsers(
            dest='convert_type',
            help='Specify the type of input file. Supported types: vcf, csv, custom.',
            required=True,
    )

    vcf2allele.setup_parser(convert_parser)

    return parser
