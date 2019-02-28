import os.path
import zipfile
import requests
import socket

'''Downloads files from google drive'''

def extract_zip(googleID, destination):
    # download zipfile in current directory, temporarily
    zipfilename = 'temp.zip'
    download(googleID, '', zipfilename)
    # extract zipfile
    with zipfile.ZipFile(zipfilename,'r') as zip_ref:
        zip_ref.extractall(destination)
    # remove zipfile
    os.remove(zipfilename)

def download(googleID, destination, filename, 
        url='https://docs.google.com/uc?export=download'):
    # path to file
    my_file = str(os.path.join(destination, filename))
    if not os.path.isfile(my_file): 
        session = requests.Session()
        assert check_internet(), \
                "HEEPS can't download input files. Check your internet connection."
        response = session.get(url, params={'id':googleID}, stream=True)
        token = get_confirm_token(response)
        if token:
            response = session.get(url, params={'id':googleID, \
                    'confirm':token}, stream=True)
        save_response_content(response, my_file)

def check_internet(host="8.8.8.8", port=53, timeout=3):
    """
    Host: 8.8.8.8 (google-public-dns-a.google.com)
    OpenPort: 53/tcp
    Service: domain (DNS/TCP)
    """
    try:
        socket.setdefaulttimeout(timeout)
        socket.socket(socket.AF_INET, socket.SOCK_STREAM).connect((host, port))
        return True
    except OSError as err:
        print('OSError [%s]: %s.'%err.args[:2])
        return False

def get_confirm_token(response):
    for key, value in response.cookies.items():
        if key.startswith('download_warning'):
            return value
    return None

def save_response_content(response, my_file, chunk_size=32768):
    with open(my_file, "wb") as f:
        for chunk in response.iter_content(chunk_size=chunk_size):
            if chunk: # filter out keep-alive new chunks
                f.write(chunk)
