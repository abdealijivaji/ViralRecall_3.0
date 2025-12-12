import logging, pyrodigal_gv, pyhmmer, pyfaidx, pandas
from pathlib import Path
#from .viralrecall import __version__
__version__ = 3.0
def setup_logger(outbase_dir : Path, name):
    
    # set up logging to file - see previous section for more details
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(message)s',
                        datefmt='%y-%m-%d %H:%M',
                        filename= outbase_dir / 'file.log',
                        filemode='w')
    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s %(message)s',   datefmt='%y-%m-%d %H:%M')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

    # Now, we can log to the root logger, or any other logger. First the root...
    logging.debug(f"Viralrecall Version: {__version__}")
    logging.debug(f'pyrodigal_gv version: {pyrodigal_gv.__version__}')
    logging.debug(f'pyhmmer version: {pyhmmer.__version__}') # type: ignore
    logging.debug(f'pyfaidx version: {pyfaidx.__version__}') # type: ignore
    logging.debug(f"Pandas version: {pandas.__version__}") 

    logger1 = logging.getLogger(name)
    
    return logger1
