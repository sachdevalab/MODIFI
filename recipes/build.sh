#!/bin/bash
# Bioconda build script for MODIFI (C++ binary + Python pipeline).
# Expects conda build env: $SRC_DIR = extracted source, $PREFIX = install prefix.

set -euo pipefail

# Compile C++ binary using the conda-provided compiler
cd "$SRC_DIR/src"
$CXX -O2 -o get_control_IPD get_control_IPD.cpp -pthread
cd "$SRC_DIR"

# Install under $PREFIX/share/modifi so that sys.path[0] = share/modifi and
# kmer_bin = share/modifi/src/get_control_IPD works without patching main.py
SHARE="$PREFIX/share/modifi"
mkdir -p "$SHARE/src"
mkdir -p "$SHARE/scripts"

# Python entry point and config loader
cp main.py load_cfg.py "$SHARE/"

# Bundled fallback for motif finding (used when pbmotifmaker is not in PATH)
if [ -d "$SRC_DIR/dependency" ]; then
  cp -r "$SRC_DIR/dependency" "$SHARE/"
fi

# Copy all pipeline scripts
cp -r "$SRC_DIR/scripts/." "$SHARE/scripts/"

# Install compiled binary
cp "$SRC_DIR/src/get_control_IPD" "$SHARE/src/"
chmod +x "$SHARE/src/get_control_IPD"

# Control database and bundled test dataset, so the layout under share/modifi
# mirrors the repo (main.py at share root, control_db/ and test/ beside it) and
# the self-locating test scripts resolve ../../main.py and ../../control_db/.
if [ -d "$SRC_DIR/control_db" ]; then
  cp -r "$SRC_DIR/control_db" "$SHARE/"
fi
if [ -d "$SRC_DIR/test" ]; then
  cp -r "$SRC_DIR/test" "$SHARE/"
fi

# Entry point wrapper so users can run "modifi" from the command line
mkdir -p "$PREFIX/bin"
cat <<'WRAPPER' > "$PREFIX/bin/modifi"
#!/bin/bash
DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
exec "$DIR/python" -B "$DIR/../share/modifi/main.py" "$@"
WRAPPER
chmod +x "$PREFIX/bin/modifi"

# "modifi-test" wrapper: one fixed command that runs the bundled HiFi test in
# the current directory (writes ./modifi_test/output/).
cp "$SRC_DIR/bin/modifi-test" "$PREFIX/bin/modifi-test"
chmod +x "$PREFIX/bin/modifi-test"
