import requests
import os

'''Downloads large mulit-cube phase screens from google drive'''

def download_from_gdrive(id, destination, filename):    
    my_file = str(destination + '/'+ filename)
    destination = destination + '/'+ filename
    if (os.path.isfile(my_file) == False): 
        print('Downloading multi-cube phase screen from the google drive.....')
        URL = "https://docs.google.com/uc?export=download"    
        session = requests.Session()    
        response = session.get(URL, params = { 'id' : id }, stream = True)
        token = get_confirm_token(response)    
        if token:
            params = { 'id' : id, 'confirm' : token }
            response = session.get(URL, params = params, stream = True)
    
        save_response_content(response, destination)    

def get_confirm_token(response):
    for key, value in response.cookies.items():
        if key.startswith('download_warning'):
            return value

    return None

def save_response_content(response, destination):
    CHUNK_SIZE = 32768

    with open(destination, "wb") as f:
        for chunk in response.iter_content(CHUNK_SIZE):
            if chunk: # filter out keep-alive new chunks
                f.write(chunk)
