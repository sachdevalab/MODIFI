import runpy
import sys
from pathlib import Path
import sysconfig


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
    Console entry point for the MODIFI pipeline.
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


if __name__ == "__main__":
    main()

