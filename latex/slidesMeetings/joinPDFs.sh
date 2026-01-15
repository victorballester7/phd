#!/bin/bash

# Color codes
RED="\e[31m"
GREEN="\e[32m"
YELLOW="\e[33m"
CYAN="\e[36m"
NC="\e[0m"

# Output PDF
output="joined_mains.pdf"

# Find all main.pdf files in directories starting with 20*, sorted
files=()
for dir in $(find . -maxdepth 1 -type d -name "20*" | sort -V); do
    if [[ -f "$dir/main.pdf" ]]; then
        files+=("$dir/main.pdf")
    fi
done

# Check if we have any PDFs to merge
if [ ${#files[@]} -eq 0 ]; then
    echo -e "${RED}Error:No main.pdf files found in directories starting with 20*.${NC}"
    exit 1
fi

# Merge using pdfunite (you can also use ghostscript or pdftk)
echo -e "${CYAN}Merging ${#files[@]} files into $output...${NC}"
pdfunite "${files[@]}" "$output"
echo -e "${GREEN}Done!${NC}"

