#!/bin/bash
# Double-click this file in Finder to launch the k-BAAM Outputs GUI.
# Requires the conda base environment (or any Python env with the dependencies).

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_PORT="$SCRIPT_DIR/pythonPort"

# Activate conda base if not already in an env
if [ -f "/opt/anaconda3/etc/profile.d/conda.sh" ]; then
    source "/opt/anaconda3/etc/profile.d/conda.sh"
    conda activate base 2>/dev/null
fi

cd "$PYTHON_PORT"
python kbaam_gui.py
