def convert_command(args):
    print("running convert commands")
    print(args)

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
    parser.add_argument(
        '--input', '-i',
        type=str,
        required=True,
        help='Path to the input file for conversion.'
    )
    parser.add_argument(
        '--output', '-o',
        type=str,
        help='Optional path for the output file. If not specified, output will be printed to stdout.'
    )

    parser.set_defaults(func=convert_command)

    return parser
