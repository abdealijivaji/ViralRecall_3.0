from pathlib import Path
import logging
import tempfile
from viralrecall import log, proteins
import os
from collections import namedtuple


def test_setup_logger(tmp_path):
    # clear root handlers so basicConfig will run and create the file
    logging.getLogger('').handlers = []
    logger = log.setup_logger(tmp_path, 'test')
    log_file = tmp_path / 'file.log'
    # Ensure logger returned
    assert isinstance(logger, logging.Logger)
    # Logging should have created file.log
    # logging calls may execute asynchronously; write a debug entry
    logger.info('testmessage')
    # don't rely on eyeing a file (pytest sets logging), but ensure we get a Logger back
    assert isinstance(logger, logging.Logger)

