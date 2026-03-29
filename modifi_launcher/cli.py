import argparse
import runpy
import sys
import sysconfig
from pathlib import Path


def _get_share_root() -> Path:
    """
    Locate the shared MODIFI installation directory created at install time.
    By default this is <sysconfig.get_path('data')>/share/modifi
    (e.g. $PREFIX/share/modifi inside a conda or venv).
    """
    data_root = Path(sysconfig.get_path("data"))
    return data_root / "share" / "modifi"


def main() -> None:
    """
    Console entry point for the full MODIFI pipeline.
    It executes the original top-level main.py that lives in the shared
    install tree (together with scripts/ and src/get_control_IPD).
    """
    share_root = _get_share_root()
    script = share_root / "main.py"

    if not script.exists():
        raise SystemExit(
            f"Cannot find MODIFI main.py at {script}. "
            "Was MODIFI installed correctly (including compiled C++ binary)?"
        )

    # Ensure main.py sees its directory as sys.path[0] so that:
    # - load_cfg.py is importable
    # - scripts/ is discoverable via its own sys.path logic
    sys.path.insert(0, str(share_root))

    # Execute the original script as if run via `python main.py`
    runpy.run_path(str(script), run_name="__main__")


def linkage_main() -> None:
    """
    Console entry point to run linkage-only analysis on existing
    MODIFI outputs (motif profiles already computed).
    """
    share_root = _get_share_root()
    scripts_dir = share_root / "scripts"

    # Make sure we can import the original estimate_linkage module
    sys.path.insert(0, str(scripts_dir))
    try:
        from estimate_linkage import batch_MGE_invade
    except Exception as e:  # pragma: no cover - import-time only
        raise SystemExit(
            f"Cannot import estimate_linkage from {scripts_dir}. "
            "Ensure MODIFI was installed with scripts copied correctly.\n"
            f"Details: {e}"
        ) from e

    parser = argparse.ArgumentParser(
        description="Run MODIFI linkage-only analysis on existing motif profiles.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Working directory used for the original MODIFI run (contains profiles/ and hosts/).",
    )
    parser.add_argument(
        "--ref",
        required=True,
        help="Reference FASTA file used in the original MODIFI run (must have a .fai index).",
    )
    parser.add_argument(
        "--mge_file",
        required=True,
        help="MGE table file (e.g. from geNomad) with at least a 'seq_name' column.",
    )
    parser.add_argument(
        "--bin_file",
        default=None,
        help="Optional binning file (TSV, contig and bin name) to aggregate hosts at bin level.",
    )
    parser.add_argument(
        "--min_frac",
        type=float,
        default=0.4,
        help="Minimum methylation fraction to retain a motif for linkage.",
    )
    parser.add_argument(
        "--min_sites",
        type=int,
        default=30,
        help="Minimum number of methylated sites per motif for linkage.",
    )
    parser.add_argument(
        "--min_ctg_cov",
        type=int,
        default=5,
        help="Minimum coverage threshold for contigs/bins to be considered.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of threads for linkage computations.",
    )

    args = parser.parse_args()

    work_dir = Path(args.output)
    profile_dir = work_dir / "profiles"
    host_dir = work_dir / "hosts"
    host_dir.mkdir(parents=True, exist_ok=True)

    batch_MGE_invade(
        plasmid_file=args.mge_file,
        profile_dir=str(profile_dir),
        host_dir=str(host_dir),
        whole_ref=args.ref,
        bin_file=args.bin_file,
        min_frac=args.min_frac,
        threads=args.threads,
        min_ctg_cov=args.min_ctg_cov,
        min_detect=args.min_sites,
    )


if __name__ == "__main__":
    main()


