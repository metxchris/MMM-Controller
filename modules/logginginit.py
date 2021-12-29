"""Initializes logging using user settings

If this module is not imported, then the default config for the logging level
is "WARNING".
"""

import logging

import settings

if settings.PRINT_SAVE_MESSAGES:
    logging.basicConfig(level="INFO")
else:
    logging.basicConfig(level="WARNING")
