#Main python script to load config file and execture running script
import sys
from pathlib import Path
from .config import load_config
from .running import run


def main() -> None:
    yaml_path = Path(__file__).parent / 'config.yaml'
    config_path = sys.argv[1] if len(sys.argv) > 1 else yaml_path   #Default to standard config file if user does not provide specific file
    
    cfg = load_config(config_path)     # parse + validate
    run(cfg)                            # orchestrate + execute
    



if __name__ == '__main__':
    main()