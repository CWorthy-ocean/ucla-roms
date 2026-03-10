#!/bin/bash

# Set the size limit in bytes (e.g., 5MB = 5242880 bytes)
MAX_SIZE=5242880

# Detect the operating system to use the correct stat syntax
if [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS
    STAT_CMD="stat -f%z"
else
    # Linux and other Unix-like systems
    STAT_CMD="stat -c%s"
fi

# Find files that exceed the size limit
large_files=$(git diff --cached --name-only --diff-filter=A | while read -r file; do
    if [ -f "$file" ]; then
        file_size=$($STAT_CMD "$file")
        if [ "$file_size" -gt "$MAX_SIZE" ]; then
            echo "$file"
        fi
    fi
done)

if [ -n "$large_files" ]; then
    echo "Error: The following files are larger than the allowed limit of $(($MAX_SIZE / 1024 / 1024))MB:"
    echo "$large_files"
    exit 1
fi

exit 0
