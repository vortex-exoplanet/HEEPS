import subprocess
import sys
import os
import heeps

def _getUserInfo():
    import getpass
    log =getpass.getuser() + '\n'
    return '----------\nUser name:\n' + log

def _getPlatformInfo():
    import platform
    log = ''.join(platform.uname()) + '\n'
    return '----------\nPlatform info:\n' + log


def _getTimeNow():
    from datetime import datetime
    dateTimeObj = datetime.now()
    log = dateTimeObj.strftime("%d-%b-%Y (%H:%M:%S.%f)") + '\n'
    return log


def _getGitInfo():
    # this_path = os.path.split(os.path.abspath(__file__))[0]
    this_path = os.path.split(heeps.__file__)[0]
    try:
        log = subprocess.check_output('git -C %s log -1'%this_path, shell=True).decode(
            sys.stdout.encoding)
    except Exception:
        log = 'not supported (very old git version)\n'
    return '----------\nGit info (HEEPS):\n' + log

def _getVipGitInfo():
    import vip_hci
    vip_path = os.path.split(vip_hci.__file__)[0]
    try:
        log = subprocess.check_output('git -C %s log -1'%vip_path, shell=True).decode(
            sys.stdout.encoding)
    except Exception:
        log = 'not available (not a git repo or very old git version)\n'
    return '----------\nGit info (VIP):\n' + log

def _getPackageVersions():
    from platform import python_version
    import proper as pr
    # import scopesim as ss
    import numpy as np
    import scipy as sp
    import matplotlib as ml
    import vip_hci as vh
    import photutils as pu
    import astropy as ap
    import skimage as si

    log = '----------\nPackage info: \n' + 'python v%s'%python_version() + '\n'
    for module in [pr, np, sp, ml, vh, pu, ap, si]:
        log += '%s v%s'%(module.__name__, module.__version__) + '\n'

    return log

def _getCUDAToolkitInfo():
    try:
        log = subprocess.check_output('conda list cudatoolkit', shell=True).decode(
            sys.stdout.encoding)
        return '----------\nCUDA (via conda) toolkit info:\n' + log
    except:
        log = subprocess.check_output('nvcc --version', shell=True).decode(
            sys.stdout.encoding)
        return '----------\nCUDA (systemwide) toolkit info:\n' + log


def _getNvidiaGPUInfo():
    log = subprocess.check_output('nvidia-smi -L', shell=True).decode(
        sys.stdout.encoding)
    return '----------\nNvidia GPU info:\n' + log


def _getNvidiaDriverInfo():
    log = subprocess.check_output('nvidia-smi', shell=True).decode(
        sys.stdout.encoding)
    lines = log.split('\n')
    return '----------\nNvidia driver info:\n' + lines[2] + '\n'


def getFingerprint():
    # fingerprint = 'METIS Simulation Environment Fingerprint:\n' + _getTimeNow() \
    #              + _getUserInfo() + _getGitInfo() + _getPlatformInfo() + _getCUDAToolkitInfo() \
    #              + _getNvidiaDriverInfo() + _getNvidiaGPUInfo()
    fingerprint = 'HEEPS Simulation Environment Fingerprint:\n' + _getTimeNow() \
                 + _getUserInfo() + _getGitInfo() + _getVipGitInfo() + _getPlatformInfo() \
                 + _getPackageVersions()
    
    return fingerprint

def writeToFile(directory, fname='fingerprint.txt'):
    fingerprint= getFingerprint()
    with open(directory + fname, "w") as file:
        file.write(fingerprint)
    return fingerprint


if __name__ == "__main__":
    print(getFingerprint())
