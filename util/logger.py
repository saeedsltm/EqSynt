import logging
import os
import shutil
import os

# Create new logger


def myLogger(logger_name, mode, level=logging.INFO):
    """
    Method to return a custom logger with the given name and level
    """
    logger = logging.getLogger(logger_name)
    logger.setLevel(level)
    format_string = ("%(asctime)s %(message)s")
    log_format = logging.Formatter(format_string)
    file_handler = logging.FileHandler("%s.log" % (logger_name), mode=mode)
    file_handler.setFormatter(log_format)
    if (logger.hasHandlers()):
        logger.handlers.clear()
    logger.addHandler(file_handler)
    return logger

# Move log file to event folder


def moveLog(eventDir):
    logFile = "%s.log" % (eventDir)
    if os.path.exists(logFile):
        shutil.copy(logFile, eventDir)
        os.remove(logFile)
