#!/bin/bash
# Bioconda build script for MODIFI (C++ binary + Python pipeline).
# Expects conda build env: $SRC_DIR = extracted source, $PREFIX = install prefix.

set -e

# Compile C++ binary (requires htslib and pbbam from conda host env)
export PKG_CONFIG_PATH="$PREFIX/lib/pkgconfig"
export LD_LIBRARY_PATH="$PREFIX/lib:$LD_LIBRARY_PATH"

cd "$SRC_DIR/src"
g++ -o get_control_IPD get_control_IPD.cpp $(pkg-config --cflags --libs htslib pbbam)
cd "$SRC_DIR"

# Install under $PREFIX/share/modifi so that sys.path[0] = share/modifi and
# kmer_bin = share/modifi/src/get_control_IPD works without patching main.py
SHARE="$PREFIX/share/modifi"
mkdir -p "$SHARE/scripts"
mkdir -p "$SHARE/src"

# Python entry point and config loader
cp main.py load_cfg.py "$SHARE/"
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
