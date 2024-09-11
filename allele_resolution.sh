#!/bin/bash

# Example usage:
# ./allele_resolution.sh one file1.tsv file2.tsv
# pypop -c  minimal.ini file1.one_fields.tsv
# pypop -c  minimal.ini file2.one_fields.tsv


# Check if at least two arguments are provided (resolution + at least one file)
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <resolution> <input_file> [<output_file>] [--prefix <prefix>]"
    echo "Valid resolutions: one, two, three"
    exit 1
fi

# Extract and validate the resolution parameter
resolution=$1
input_file=$2
prefix='\*'

# Process optional arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --prefix)
            prefix=$2
            shift 2
            ;;
        *)
            shift
            ;;
    esac
done



case $resolution in
    one)
        suffix=".1field.tsv"
        ;;
    two)
        suffix=".2fields.tsv"
        ;;
    three)
        suffix=".3fields.tsv"
        ;;
    *)
        echo "Invalid resolution: $resolution. Valid options are: one, two, three."
        exit 1
        ;;
esac

output_file=${3:-$(basename "$input_file")$suffix}

# Check if the input file exists and is readable
if [ ! -f "$input_file" ] || [ ! -r "$input_file" ]; then
    echo "File $input_file does not exist or is not readable."
    exit 1
fi

# Check if the input file is a TSV file
if [[ $input_file != *.tsv ]]; then
    echo "File $input_file is not a TSV file. The file name should end with .tsv"
    exit 1
fi

## Filter by resolution

# Define regex patterns
one_field="($prefix[0-9]{2})(\t|$)"
two_field="($prefix[0-9]{2})(:[0-9]{2})(\t|$)"
three_field='(:[0-9]{2}){2,}G?'

# Apply the corresponding sed command based on the resolution
case $resolution in
    one)
        sed -E -e ":a; s/$two_field/\1\3/g; ta" "$input_file" | sed -E "s/$three_field//g" > "$output_file"
        ;;
    two)
        sed -E ":a; s/$one_field/\1:01\2/g; ta" "$input_file" | sed -E "s/$three_field/\1/g" > "$output_file"
        ;;
    three)
        sed -E ":a; s/$one_field/\1:01:01\2/g; ta" "$input_file" | sed -E ":a; s/$two_field/\1\2:01\3/g; ta" > "$output_file"
        ;;
esac