import os
import glob
import shutil

# https://stackoverflow.com/a/5423147
_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__),"../","../"))
def get_gendir():
    return os.path.join(_ROOT, 'generated')
def get_inpdir():
    return os.path.join(_ROOT, 'src','thermotools','data')

def empty_dir(outdir):
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir)
