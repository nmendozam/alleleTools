import vcf
import sys

if __name__ == '__main__':
    file = sys.argv[1]
    out = sys.argv[2]
    vcf_reader = vcf.Reader(filename=file)
    vcf_writer = vcf.Writer(open(out, 'w'), vcf_reader)
    for record in vcf_reader:
        vcf_writer.write_record(record)
