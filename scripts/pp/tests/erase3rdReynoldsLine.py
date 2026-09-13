#!/usr/bin/env python3
"""
Remove the first line of each group of three repeated Reynolds numbers.
Keeps only the last two lines of each group of three.
"""

def process_file(input_file, output_file=None):
    """
    Process the file: when three consecutive lines have the same first column (Re),
    remove the first line of that group.
    """
    if output_file is None:
        output_file = input_file
    
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # Skip header lines (lines starting with '#')
    header_lines = []
    data_lines = []
    for line in lines:
        if line.startswith('#'):
            header_lines.append(line)
        else:
            data_lines.append(line)
    
    # Process data lines
    filtered_lines = []
    i = 0
    while i < len(data_lines):
        # Get current Reynolds number
        current_re = data_lines[i].split()[0]
        
        # Check if there are at least 3 lines with the same Reynolds number
        if i + 2 < len(data_lines):
            next_re = data_lines[i+1].split()[0]
            third_re = data_lines[i+2].split()[0]
            
            if current_re == next_re == third_re:
                # Skip the first line (i), keep lines i+1 and i+2
                filtered_lines.append(data_lines[i+1])
                filtered_lines.append(data_lines[i+2])
                i += 3
                continue
        
        # If not a group of 3, keep the line
        filtered_lines.append(data_lines[i])
        i += 1
    
    # Write output
    with open(output_file, 'w') as f:
        f.writelines(header_lines)
        f.writelines(filtered_lines)
    
    print(f"Processed {len(data_lines)} lines -> {len(filtered_lines)} lines")
    print(f"Removed {len(data_lines) - len(filtered_lines)} lines")
    print(f"Output written to: {output_file}")
    
    return output_file


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python remove_first_re.py <input_file> [output_file]")
        print("Example: python remove_first_re.py data.txt")
        print("         python remove_first_re.py data.txt output.txt")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    process_file(input_file, output_file)
