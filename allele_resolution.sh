#!/bin/bash

# Example usage:
# ./allele_resolution.sh one file1.tsv file2.tsv
# pypop -c  minimal.ini file1.one_fields.tsv
# pypop -c  minimal.ini file2.one_fields.tsv

# Define regex patterns
one_field='(\t[0-9]{2})(\t|$)'
two_field='(\t[0-9]{2})(:[0-9]{2})(\t|$)'
three_field='(:[0-9]{2}){2,}G?'

# Check if at least two arguments are provided (resolution + at least one file)
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <resolution> <file1> [<file2> ...]"
    echo "Valid resolutions: one, two, three"
    exit 1
fi

# Extract and validate the resolution parameter
resolution=$1
case $resolution in
    one|two|three)
        ;;
    *)
        echo "Invalid resolution: $resolution. Valid options are: one, two, three."
        exit 1
        ;;
esac

# Iterate over the list of files
for file in "${@:2}"; do
    if [ ! -f "$file" ] || [ ! -r "$file" ]; then
        echo "File $file does not exist or is not readable."
        continue
    fi
    # check file format (tsv)
    if [[ $file != *.tsv ]]; then
        echo "File $file is not a TSV file. The file name should end with .tsv"
        continue
    fi

    base=$(basename "$file")

    # Apply the corresponding sed command based on the resolution
    case $resolution in
        one)
            sed -E -e ":a; s/$two_field/\1\3/g; ta" "$file" | sed -E "s/$three_field//g" > "${base}.1field.tsv"
            ;;
        two)
            sed -E ":a; s/$one_field/\1:01\2/g; ta" "$file" | sed -E "s/$three_field/\1/g" > "${base}.2fields.tsv"
            ;;
        three)
            sed -E ":a; s/$one_field/\1:01:01\2/g; ta" "$file" | sed -E ":a; s/$two_field/\1\2:01\3/g; ta" > "${base}.3fields.tsv"
            ;;
    esac

    echo "Processed $file into ${file}.${resolution}_fields.tsv"
done
