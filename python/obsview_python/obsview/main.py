#Main python script to load config file and execture running script
import sys
from pathlib import Path

from.cli import parse_args
from .config import load_config
from .reference import load_reference
from .running import run


def main() -> None:
    args = parse_args()
    
    cfg = load_config(args.config)
    ref = load_reference(args.reference)
    
    run(cfg, ref)
    



if __name__ == '__main__':
    main()