import tarfile
import zipfile
import requests
import os
import hashlib

# Calculate the checksum of a file using the BLAKE2b algorithm
def checksum(filename:str):
    # Adapted from https://stackoverflow.com/a/1131238
    with open(filename, "rb") as hdl:
        file_hash = hashlib.blake2b()
        while chunk := hdl.read(8192):
            file_hash.update(chunk)
    return file_hash.hexdigest()

# Write checksum for a file
def writesum(filename:str):
    chkpath = filename+".chk"
    print("    checksum: %s"%chkpath)
    with open(chkpath, "w") as hdl:
       hdl.write("%s \n"%checksum(filename))

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

