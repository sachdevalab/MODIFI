#!/bin/bash
# Bioconda build script for MODIFI (C++ binary + Python pipeline).
# Expects conda build env: $SRC_DIR = extracted source, $PREFIX = install prefix.

set -e

# Compile C++ binary (only needs standard C++ and pthreads)
cd "$SRC_DIR/src"
g++ -O2 -o get_control_IPD get_control_IPD.cpp -pthread
cd "$SRC_DIR"

# Install under $PREFIX/share/modifi so that sys.path[0] = share/modifi and
# kmer_bin = share/modifi/src/get_control_IPD works without patching main.py
SHARE="$PREFIX/share/modifi"
mkdir -p "$SHARE/scripts"
mkdir -p "$SHARE/src"

# Python entry point and config loader
cp main.py load_cfg.py "$SHARE/"
if [ -d "$SRC_DIR/dependency" ]; then
  cp -r "$SRC_DIR/dependency" "$SHARE/"
fi
cp -r scripts/*.py "$SHARE/scripts/"
cp src/get_control_IPD "$SHARE/src/"
chmod +x "$SHARE/src/get_control_IPD"

# Entry point script so users can run "modifi" from the command line
mkdir -p "$PREFIX/bin"
# Use dir of this script so "modifi" works even if CONDA_PREFIX is not set
cat <<'WRAPPER' > "$PREFIX/bin/modifi"
#!/bin/bash
DIR=$(cd "$(dirname "$0")" && pwd)
exec python -B "$DIR/../share/modifi/main.py" "$@"
WRAPPER
chmod +x "$PREFIX/bin/modifi"
