import logging

import settings

if settings.PRINT_SAVE_MESSAGES:
    logging.basicConfig(level="INFO")
else:
    logging.basicConfig(level="WARNING")
