#Function to parse(break down) command line arguments

import argparse
from pathlib import Path

def parse_args() -> argparse.Namespace:
    base_dir = Path(__file__).parent
    
    parser = argparse.ArgumentParser(description="Run script with config and reference files.")
    parser.add_argument(
        "-c", "--config",
        type=Path,
        default=base_dir / "config.yaml",
        help="Path to the config YAML file (default: config.yaml)"
    )
    parser.add_argument(
        "-r", "--reference",
        type=Path,
        default=base_dir / "reference.yaml",
        help="Path to the reference YAML file (default: reference.yaml)"
    )
    
    return parser.parse_args()