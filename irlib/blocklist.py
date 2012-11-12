"""
This file keeps track of which radar lines are excluded from analyses because
of missing or corrupted data.
"""

import os.path as path

def blocklist(survey_path, line):
    """ Takes a survey filename and a line number, and compares to an index to
    determine whether or not the line should be blocked. Returns a boolean.
    """
    survey_name = path.basename(survey_name)
    
    if survey_name == "glacier1_08_utm.h5":
        if line <= 9:
            block = True
    
    else:
        block = False

    return block

