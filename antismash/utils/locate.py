'''Helper function for locating files'''
import logging
import os
from os import path
import sys


def get_full_path(current_file, file_to_add):
    "Get the full path of file_to_add in the same directory as current_file"
    return path.join(path.dirname(path.abspath(current_file)), file_to_add)


def locate_executable(name):
    "Find an executable in the path and return the full path"
    # In windows, executables tend to end on .exe
    if sys.platform == 'win32':
        name += ".exe"
    file_path, _ = os.path.split(name)
    if file_path != "":
        if path.isfile(name) and os.access(name, os.X_OK):
            logging.debug("Found executable %r", name)
            return name
    for env_path in os.environ["PATH"].split(os.pathsep):
        full_name = path.join(env_path, name)
        if path.isfile(full_name) and os.access(full_name, os.X_OK):
            logging.debug("Found executable %r", full_name)
            return full_name

    return None


def locate_file(name):
    "Find a file and return the full path"
    file_path, _ = os.path.split(name)
    if file_path != "":
        if path.isfile(name) and os.access(name, os.R_OK):
            logging.debug("Found file %r", name)
            return name
    return None
