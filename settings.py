# MMM Driver path (USE FORWARDSLASHES)
MMM_DRIVER_PATH = 'C:/cygwin64/home/metxc/mmm/wrapper/mmm.exe'
# MMM_DRIVER_PATH = 'C:/cygwin64/home/metxc/mmm/test/mmm.exe'

# Make Profile PDFS when running scans
MAKE_PROFILE_PDFS = True

# Automatically open PDFs after they are merged
AUTO_OPEN_PDFS = True

# Print messages for all saved files
PRINT_SAVE_MESSAGES = False

# Print non-error responses from MMM
PRINT_MMM_RESPONSE = False

# Reduce file sizes by not saving additional variables
SAVE_ADDITIONAL_VARIABLES = True

# Adjust the initial scan number to group different types of scans
STARTING_SCAN_NUMBER = 1

# Method of interpolation using scipy.interpolate.interp1d(kind=INTERPOLATION_METHOD)
#   Consider using additional smoothing with linear interpolation, and less
#   smoothing with cubic interpolation.  Only used when GRADIENT_METHOD='interpolate'
#   kind: 'slinear', 'quadratic', 'cubic'
INTERPOLATION_METHOD = 'quadratic'

# Method to take gradients:
#   'akima':       This is the same gradient method used by TRANSP for MMM calls.
#                  Unfortunately, this method can be slow in Python.
#   'traditional': Gradients taken using a simple nearest neighbor method.
#                  This is faster, but less accurate when not many input points are available.
#   'ptsolver':    A worse implementation of the nearest neighbor method.
#   'interpolate': Gradients taken using interpolation.  This is more accurate when not using
#                  many input points, but the endpoint values are more dependent on the
#                  interpolation method.
GRADIENT_METHOD = 'akima'
SOLVER_GRADIENT_METHOD = 'ptsolver'

# Temp EPM switch
USE_EPM = True

# 'old', '#90', '#105', '#107', '#111', '#113', '#114', '#117', '#123', '#129'
MMM_HEADER_VERSION = '#123'
