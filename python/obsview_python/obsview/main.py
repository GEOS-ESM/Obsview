#Main python script to load config file and execture running script
import sys
from pathlib import Path
from .config import load_config
from .reference import load_reference
from .running import run


def main() -> None:
    #Load config
    yaml_path = Path(__file__).parent / 'config.yaml'
    config_path = sys.argv[1] if len(sys.argv) > 1 else yaml_path   #Default to standard config file if user does not provide specific file
    cfg = load_config(config_path)     

    #Load reference
    ref_path = Path(__file__).parent / 'reference.yaml'
    ref = load_reference(ref_path)

    run(cfg, ref)                            
    



if __name__ == '__main__':
    main()