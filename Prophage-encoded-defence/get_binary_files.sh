#!/bin/bash

# Set the destination directory for files with 200+ lines
DEST_DIR="large_binary_files"

# Create the destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

# Find all .binary files recursively and process them
find . -type f -name "*.binary" | while read -r file; do
    # Count the number of non-empty lines in the file
    line_count=$(grep -cve '^\s*$' "$file")
    
    # Check if the file has 200 or more lines
    if [[ "$line_count" -ge 200 ]]; then
        echo "File '$file' has $line_count lines. Moving to $DEST_DIR."
        mv "$file" "$DEST_DIR"
    else
        echo "File '$file' has $line_count lines. Skipping."
    fi
done

echo "Processing complete. Check the '$DEST_DIR' directory for results."
