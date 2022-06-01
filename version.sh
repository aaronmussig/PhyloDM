#!/bin/bash
set -e

# Set variables
VERSION=$1
RELEASE_DATE=$(date +"%Y-%m-%d")

# Update Cargo version
sed -i "s/^version = .*$/version = \"$VERSION\"/" Cargo.toml

# Update citation version and date
sed -i "s/^version: .*$/version: $VERSION/" CITATION.cff
sed -i "s/^date-released: .*$/date-released: $RELEASE_DATE/" CITATION.cff

