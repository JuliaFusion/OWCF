import sys
import os
PYTHON_FOLDERPATH = os.path.dirname(sys.executable)
with open('PYTHON_FOLDERPATH.txt', 'w') as f:
    f.write(PYTHON_FOLDERPATH)