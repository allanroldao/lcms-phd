#!/bin/bash

# Directory containing the files
directory="/Users/allanpradelliroldao/Documents/Covalent_Fragments/raw_data/both_equipments/pools"

# Output file
output_file="merged_data_within_interval.txt"

# Define the interval of x-axis values
x_min=29300
x_max=30200

# Clear the output file if it exists
> "${output_file}"

# Loop through all files in the directory
for i in {0..160}; do
    file_name="ncSH2003_pool${i}.txt"
    file_path="${directory}/${file_name}"
    
    # Check if the file exists
    if [ -f "$file_path" ]; then
        echo "Processing: ${file_name}"

        echo "Sample: ${file_name}" >> "${output_file}"
        
        # Determine the delimiter based on file name
        if [ $i -eq 0 ] || [ $i -ge 58 ]; then
            delimiter=","
            skip_lines=1
        else
            delimiter="\t"
            skip_lines=0
        fi
        
        # Filter data within the x-axis interval, using only the integer part of the first column
        if [ "$delimiter" == "," ]; then
            # For comma-separated files with headers
            awk -F"$delimiter" -v OFS="\t" -v x_min="$x_min" -v x_max="$x_max" -v skip_lines="$skip_lines" \
                'NR > skip_lines {
                    gsub(/"/, "", $1);    # Remove quotes from the first column
                    gsub(/"/, "", $2);    # Remove quotes from the second column
                    mass_int=int($1);     # Convert first column to integer (ignore decimals)
                    if (mass_int >= x_min && mass_int <= x_max) {
                        print mass_int, $2;  # Print integer mass and second column
                    }
                }' "${file_path}" >> "${output_file}"
        else
            # For tab-separated files without headers
            awk -F"$delimiter" -v OFS="\t" -v x_min="$x_min" -v x_max="$x_max" -v skip_lines="$skip_lines" \
                'NR > skip_lines {
                    mass_int=int($1);     # Convert first column to integer (ignore decimals)
                    if (mass_int >= x_min && mass_int <= x_max) {
                        print mass_int, $2;  # Print integer mass and second column
                    }
                }' "${file_path}" >> "${output_file}"
        fi
        
        echo "" >> "${output_file}"  # Add a blank line after each sample
    else
        echo "File ${file_name} does not exist."
    fi
done
