import os
import sys
import yaml
import shutil
from pathlib import Path

# Directory containing this module (repo root or pip/share/modifi layout). Do not use
# sys.path[0] for bundled assets — it is not guaranteed to match after imports.
_MODIFI_ROOT = Path(__file__).resolve().parent


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
    motif_maker_bin = None
    pbmm2_bin = None
    pbindex_bin = None
    
    # 1️⃣ Try config.yaml first (alongside main.py / this module)
    config_file = _MODIFI_ROOT / "smrt.config.yaml"
    print (f"Config file: {config_file}")
    if config_file.exists():
        with open(config_file, "r") as cf:
            config = yaml.safe_load(cf)
        smrt_bin = Path(config.get("smrtlink_bin", ""))
        
        if smrt_bin.exists():
            # Try to find tools in SMRT Link directory
            motif_maker_candidate = smrt_bin / "pbmotifmaker"
            if not motif_maker_candidate.exists():
                motif_maker_candidate = smrt_bin / "motifmaker"  # fallback name
            if motif_maker_candidate.exists():
                motif_maker_bin = str(motif_maker_candidate)
            
            pbmm2_candidate = smrt_bin / "pbmm2"
            if pbmm2_candidate.exists():
                pbmm2_bin = str(pbmm2_candidate)
            
            pbindex_candidate = smrt_bin / "pbindex"
            if pbindex_candidate.exists():
                pbindex_bin = str(pbindex_candidate)
            
            if all([motif_maker_bin, pbmm2_bin, pbindex_bin]):
                print(f"✅ Using binaries from config file: {config_file}")
                print(f"  pbmotifmaker: {motif_maker_bin}")
                print(f"  pbmm2: {pbmm2_bin}")
                print(f"  pbindex: {pbindex_bin}")
                return motif_maker_bin, pbmm2_bin, pbindex_bin
        else:
            print(f"⚠️  SMRT Link path in config.yaml does not exist: {smrt_bin}")
    
    # 2️⃣ Fall back to system PATH
    print("🔍 Config file not found or incomplete, checking system PATH...")
    if not motif_maker_bin:
        motif_maker_bin = find_tool("pbmotifmaker", alt_names=["motifmaker"])
    if not pbmm2_bin:
        pbmm2_bin = find_tool("pbmm2")
    if not pbindex_bin:
        pbindex_bin = find_tool("pbindex")

    # 3️⃣ If still not all found, use MultiMotifMaker as fallback for motif maker
    if not motif_maker_bin:
        print("⚠️  pbmotifmaker not found, using MultiMotifMaker.jar instead")
        motif_maker_bin = str(_MODIFI_ROOT / "dependency" / "MultiMotifMaker.jar")

    if all([motif_maker_bin, pbmm2_bin, pbindex_bin]):
        print("✅ Using binaries found in system PATH.")
        print(f"  pbmotifmaker: {motif_maker_bin}")
        print(f"  pbmm2: {pbmm2_bin}")
        print(f"  pbindex: {pbindex_bin}")
        return motif_maker_bin, pbmm2_bin, pbindex_bin
    else:
        print("❌ Required tools not found. Please ensure motifmaker, pbmm2 and pbindex are installed and in your PATH.")

    # 4️⃣ Validate existence of final paths
    for tool_name, tool_path in [("pbmotifmaker", motif_maker_bin), ("pbmm2", pbmm2_bin), ("pbindex", pbindex_bin)]:
        if not Path(tool_path).exists():
            raise FileNotFoundError(f"❌ Required tool not found: {tool_name} at {tool_path}")

    return motif_maker_bin, pbmm2_bin, pbindex_bin


# Example usage
if __name__ == "__main__":
    motif_maker_bin, pbmm2_bin, pbindex_bin = load_binaries()

