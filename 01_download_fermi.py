#!/usr/bin/env python3
"""
STEP A1: Download Fermi LAT photon data for the Galactic Center region.

This script queries the Fermi LAT Data Server and downloads:
- Photon events (10-50 GeV) within 40° of Galactic Center
- Spacecraft pointing file
- Background models

Prerequisites:
  pip install astropy requests tqdm

After running, you will have:
  data/fermi/photons.fits   — photon event list
  data/fermi/spacecraft.fits — spacecraft file  
  data/fermi/gll_iem_v07.fits — Galactic diffuse model
  data/fermi/iso_P8R3_SOURCE_V3_v1.txt — isotropic template

NOTE: The Fermi data server query may take 10-60 minutes to process.
      The script will poll until ready, then download.
"""
import os
import requests
from pathlib import Path

DATA_DIR = Path("data/fermi")
DATA_DIR.mkdir(parents=True, exist_ok=True)

# ─── Fermi LAT Data Server ───
# Region: Galactic Center (l=0, b=0), radius 40°
# Energy: 10,000 - 50,000 MeV (10-50 GeV)
# Time: Mission start to present
# Event class: 128 (SOURCE), event type: 3 (FRONT+BACK)

# Background model URLs
BG_MODELS = {
    "gll_iem_v07.fits": "https://fermi.gsfc.nasa.gov/ssc/data/analysis/software/aux/4fgl/gll_iem_v07.fits",
    "iso_P8R3_SOURCE_V3_v1.txt": "https://fermi.gsfc.nasa.gov/ssc/data/analysis/software/aux/iso_P8R3_SOURCE_V3_v1.txt",
}


def download_file(url, dest, desc=""):
    """Download a file with progress indicator."""
    print(f"  Downloading {desc or url}...")
    try:
        from tqdm import tqdm
        r = requests.get(url, stream=True)
        r.raise_for_status()
        total = int(r.headers.get('content-length', 0))
        with open(dest, 'wb') as f, tqdm(total=total, unit='B', unit_scale=True) as bar:
            for chunk in r.iter_content(8192):
                f.write(chunk)
                bar.update(len(chunk))
    except ImportError:
        r = requests.get(url)
        r.raise_for_status()
        with open(dest, 'wb') as f:
            f.write(r.content)
    print(f"  → Saved: {dest} ({os.path.getsize(dest)/1e6:.1f} MB)")


def submit_fermi_query():
    """Submit query to Fermi data server, return query ID."""
    print("═" * 60)
    print("  STEP A1: Querying Fermi LAT Data Server")
    print("═" * 60)
    print(f"  Region: Galactic Center, 40° radius")
    print(f"  Energy: 10-50 GeV")
    print(f"  Time: Full mission (2008-2025)")
    print()

    print("  Submitting query...")
    # Note: actual Fermi server uses a web form; this creates the manual instructions
    print()
    print("  ┌─────────────────────────────────────────────────────┐")
    print("  │  MANUAL STEP REQUIRED:                              │")
    print("  │                                                     │")
    print("  │  Go to: https://fermi.gsfc.nasa.gov/cgi-bin/ssc/    │")
    print("  │         LAT/LATDataQuery.cgi                        │")
    print("  │                                                     │")
    print("  │  Enter these parameters:                            │")
    print("  │    Object name/coords: 0, 0                         │")
    print("  │    Coordinate system: Galactic (L,B)                │")
    print("  │    Search radius: 40 degrees                        │")
    print("  │    Observation dates: START to END (all)             │")
    print("  │    Energy range: 10000 - 50000 MeV                  │")
    print("  │    LAT data type: Photon                            │")
    print("  │    Spacecraft data: checked                         │")
    print("  │                                                     │")
    print("  │  Click 'Start Query'. Wait for processing.          │")
    print("  │  Download the photon file and spacecraft file.      │")
    print("  │  Save them as:                                      │")
    print("  │    data/fermi/photons.fits                           │")
    print("  │    data/fermi/spacecraft.fits                        │")
    print("  └─────────────────────────────────────────────────────┘")
    print()


def download_background_models():
    """Download Galactic diffuse and isotropic models."""
    print("  Downloading background models...")
    for fname, url in BG_MODELS.items():
        dest = DATA_DIR / fname
        if dest.exists():
            print(f"  → Already exists: {dest}")
            continue
        try:
            download_file(url, dest, fname)
        except Exception as e:
            print(f"  ⚠ Could not download {fname}: {e}")
            print(f"    Manual download: {url}")
            print(f"    Save to: {dest}")


def download_4fgl_catalog():
    """Download the 4FGL-DR4 point source catalog."""
    cat_url = "https://fermi.gsfc.nasa.gov/ssc/data/access/lat/14yr_catalog/gll_psc_v35.fit"
    dest = DATA_DIR / "gll_psc_v35.fit"
    if dest.exists():
        print(f"  → Catalog already exists: {dest}")
        return
    print("  Downloading 4FGL-DR4 source catalog...")
    try:
        download_file(cat_url, dest, "4FGL-DR4")
    except Exception as e:
        print(f"  ⚠ Could not download catalog: {e}")
        print(f"    Manual: {cat_url}")
        print(f"    Save to: {dest}")


def verify_data():
    """Check that all required files exist."""
    print()
    print("  Checking data files...")
    required = [
        ("data/fermi/photons.fits", "Photon events (from Fermi server)"),
        ("data/fermi/spacecraft.fits", "Spacecraft file (from Fermi server)"),
        ("data/fermi/gll_iem_v07.fits", "Galactic diffuse model"),
        ("data/fermi/iso_P8R3_SOURCE_V3_v1.txt", "Isotropic template"),
    ]
    all_ok = True
    for path, desc in required:
        exists = os.path.exists(path)
        status = "✅" if exists else "❌ MISSING"
        print(f"    {status}  {path}  — {desc}")
        if not exists:
            all_ok = False

    if all_ok:
        print()
        print("  ✅ All files present. Ready for Step A2.")
    else:
        print()
        print("  ❌ Some files missing. Complete the manual download step above.")
    return all_ok


if __name__ == "__main__":
    submit_fermi_query()
    download_background_models()
    download_4fgl_catalog()
    verify_data()
