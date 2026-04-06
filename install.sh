#!/bin/bash
# Install script: adds pipeline bin/ directory to user's PATH

set -euo pipefail

# Get the absolute path of the script's directory
SCRIPT_DIR="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"
BIN_DIR="$SCRIPT_DIR/bin"

echo "Installing DADA2 Pipeline..."
echo "Bin directory: $BIN_DIR"

# Check if PATH already contains this directory
if grep -q "$BIN_DIR" ~/.bashrc 2>/dev/null; then
    echo "✓ Already in PATH"
else
    echo "Adding $BIN_DIR to ~/.bashrc"
    cat >> ~/.bashrc << EOF

# Added by DADA2 Pipeline installer
export PATH="$BIN_DIR:\$PATH"
EOF
    echo "✓ Added to PATH"
fi

# Make sure the script is executable
chmod +x "$BIN_DIR/run_dada2processing"
echo "✓ Made run_dada2processing executable"

echo ""
echo "=========================================="
echo "Installation complete!"
echo "=========================================="
echo ""
echo "To activate the changes, run:"
echo "  source ~/.bashrc"
echo ""
echo "Then you can use:"
echo "  run_dada2processing --help"
echo ""
