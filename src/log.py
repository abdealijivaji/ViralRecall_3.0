import logging, pyrodigal_gv, pyhmmer
from pathlib import Path
__version__ = 3.0

def setup_logger(outbase_dir):
    # import logging

    # set up logging to file - see previous section for more details
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(message)s',
                        datefmt='%y-%m-%d %H:%M',
                        filename= outbase_dir / 'file.log',
                        filemode='w')
    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(asctime)s %(message)s',   datefmt='%y-%m-%d %H:%M')
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)

    # Now, we can log to the root logger, or any other logger. First the root...
    logging.debug(f"Viralrecall Version: {__version__}")
    logging.debug(f'pyrodigal_gv version: {pyrodigal_gv.__version__}')
    logging.debug(f'pyhmmer version: {pyhmmer.__version__}') # type: ignore

    # Now, define a couple of other loggers which might represent areas in your
    # application:

    logger1 = logging.getLogger('myapp.area1')
    # logger2 = logging.getLogger('myapp.area2')

    # logger1.debug('Quick zephyrs blow, vexing daft Jim.')
    # logger1.info('How quickly daft jumping zebras vex.')
    # logger2.warning('Jail zesty vixen who grabbed pay from quack.')
    # logger2.error('The five boxing wizards jump quickly.')
    return logger1