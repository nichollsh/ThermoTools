import tarfile
import zipfile
import requests
import os

# https://stackoverflow.com/a/53101953
def download(url, fpath):
    print("Download file from "+url)
    r = requests.get(url, stream=True)
    with open(fpath, 'wb') as f:
        for chunk in r.iter_content(chunk_size=1024):
            if chunk: # filter out keep-alive new chunks
                f.write(chunk)
    print("    Done")

def untar(fpath, dpath):
    print("Untarring "+fpath)
    with tarfile.open(fpath)  as hdl:
        hdl.extractall(dpath) 
    print("    Done")

def unzip(fpath, dpath):
    print("Unzipping "+fpath)
    with zipfile.ZipFile(fpath,"r") as hdl:
        hdl.extractall(dpath)
    print("    Done")

def makezip(zpath, files):
    print("Zipping to "+zpath)
    if os.path.exists(zpath):
        print("    removing old zip file")
        os.remove(zpath)
    with zipfile.ZipFile(zpath, 'w', zipfile.ZIP_DEFLATED) as hdl:
        for f in files:
            hdl.write(f, arcname=os.path.basename(f))
    print("Done")

