import os
import sys
import yaml
import shutil
from pathlib import Path

def find_tool(name, alt_names=None):
    """Try to find a tool in PATH, or from alternative names."""
    alt_names = alt_names or []
    # Check main name first
    path = shutil.which(name)
    if path:
        return path
    # Then check alternatives
    for alt in alt_names:
        path = shutil.which(alt)
        if path:
            return path
    return None

def load_binaries():
    # 1️⃣ Try system PATH first
    motif_maker_bin = find_tool("pbmotifmaker", alt_names=["motifmaker"])
    pbmm2_bin = find_tool("pbmm2")
    pbindex_bin = find_tool("pbindex")

    if all([motif_maker_bin, pbmm2_bin, pbindex_bin]):
        print("✅ Using binaries found in system PATH.")
        return motif_maker_bin, pbmm2_bin, pbindex_bin

    # 2️⃣ If not all found, look for config.yaml
    config_file = Path(sys.path[0]) / "config.yaml"
    if config_file.exists():
        with open(config_file, "r") as cf:
            config = yaml.safe_load(cf)
        smrt_bin = Path(config.get("smrtlink_bin", ""))
        if not smrt_bin.exists():
            raise FileNotFoundError(f"SMRT Link path not found: {smrt_bin}")

        motif_maker_bin = motif_maker_bin or str(smrt_bin / "pbmotifmaker")
        if not Path(motif_maker_bin).exists():
            motif_maker_bin = str(smrt_bin / "motifmaker")  # fallback name
        pbmm2_bin = pbmm2_bin or str(smrt_bin / "pbmm2")
        pbindex_bin = pbindex_bin or str(smrt_bin / "pbindex")
        print(f"✅ Loaded binaries from {config_file}")
    else:
        raise FileNotFoundError(
            "SMRT Link tools not found in PATH and no config.yaml detected. "
            "Please set PATH or create a config.yaml file with:\n"
            "smrtlink_bin: /path/to/smrtlink/private/bin"
        )

    # 3️⃣ Validate existence
    for tool in [motif_maker_bin, pbmm2_bin, pbindex_bin]:
        if not Path(tool).exists():
            raise FileNotFoundError(f"❌ Required tool not found: {tool}")

    print(
        f"Loaded binaries:\n"
        f"  pbmotifmaker: {motif_maker_bin}\n"
        f"  pbmm2: {pbmm2_bin}\n"
        f"  pbindex: {pbindex_bin}"
    )

    return motif_maker_bin, pbmm2_bin, pbindex_bin


# Example usage
if __name__ == "__main__":
    motif_maker_bin, pbmm2_bin, pbindex_bin = load_binaries()

