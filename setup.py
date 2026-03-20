from __future__ import annotations

import os
import shutil
import subprocess
import sysconfig
from pathlib import Path

from setuptools import find_packages, setup
from setuptools.command.install import install


class InstallWithBinary(install):
    """
    Custom install command that:
    1) Compiles src/get_control_IPD (standard C++ with pthreads)
    2) Installs the original Python pipeline layout into <prefix>/share/modifi
       so that main.py and scripts/ remain unchanged.
    """

    def run(self) -> None:
        super().run()
        self._build_and_install_modifi_tree()

    def _build_and_install_modifi_tree(self) -> None:
        root = Path(__file__).resolve().parent
        src_dir = root / "src"
        cpp = src_dir / "get_control_IPD.cpp"

        if not cpp.exists():
            print("WARNING: src/get_control_IPD.cpp not found; skipping C++ build.")
            return

        exe = src_dir / "get_control_IPD"

        # Compile the C++ binary by delegating to src/install.sh,
        # so we keep a single source of truth for compiler flags.
        try:
            subprocess.check_call(["bash", "install.sh"], cwd=str(src_dir))
        except Exception as e:  # pragma: no cover - build-time only
            print(
                "WARNING: failed to build get_control_IPD via src/install.sh; "
                "ensure g++ is available.",
                flush=True,
            )
            print(f"Details: {e}", flush=True)
            return

        # Install into <data>/share/modifi to match the Bioconda layout
        data_root = Path(sysconfig.get_path("data"))
        share_root = data_root / "share" / "modifi"
        scripts_dst = share_root / "scripts"
        src_dst = share_root / "src"

        share_root.mkdir(parents=True, exist_ok=True)
        scripts_dst.mkdir(parents=True, exist_ok=True)
        src_dst.mkdir(parents=True, exist_ok=True)

        # Copy core Python entry files
        shutil.copy2(root / "main.py", share_root / "main.py")
        shutil.copy2(root / "load_cfg.py", share_root / "load_cfg.py")

        # Copy all helper scripts
        scripts_src = root / "scripts"
        if scripts_src.is_dir():
            for py in scripts_src.glob("*.py"):
                shutil.copy2(py, scripts_dst / py.name)

        # Copy compiled binary
        shutil.copy2(exe, src_dst / "get_control_IPD")
        os.chmod(src_dst / "get_control_IPD", 0o755)

        print(f"Installed MODIFI pipeline tree under {share_root}", flush=True)


setup(
    name="modifi",
    version="0.0.0",
    description="DNA modification detection and MGE-host linkage inference",
    packages=find_packages(include=["modifi_launcher", "modifi_launcher.*"]),
    include_package_data=True,
    cmdclass={"install": InstallWithBinary},
    entry_points={
        "console_scripts": [
            "modifi=modifi_launcher.cli:main",
            "MODIFI=modifi_launcher.cli:main",
            "modifi-linkage=modifi_launcher.cli:linkage_main",
        ]
    },
    python_requires=">=3.8",
)

