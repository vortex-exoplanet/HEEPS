import requests
import os.path

'''Downloads files from google drive'''

def download_from_gdrive(id, destination, filename, 
        url='https://docs.google.com/uc?export=download'):
    # path to file
    my_file = str(os.path.join(destination, filename))
    if not os.path.isfile(my_file): 
        print('Downloading "%s" from google drive' %filename)
        session = requests.Session()
        response = session.get(url, params={'id':id}, stream=True)
        token = get_confirm_token(response)
        if token:
            response = session.get(url, params={'id':id, 'confirm':token}, stream=True)
        save_response_content(response, my_file)

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
