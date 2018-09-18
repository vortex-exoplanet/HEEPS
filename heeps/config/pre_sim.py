import os
import proper
import zipfile

from .download_from_gdrive import download_from_gdrive

def pre_sim(conf):
	# Creates required directories for simulation 
	os.makedirs(conf['OUT_DIR'] , exist_ok=True)
	os.makedirs(conf['TMP_DIR'], exist_ok=True)
	os.makedirs(conf['INPUT_DIR'], exist_ok=True)
	proper.print_it = False   # To suppress the printing of intermediate steps by PROPER routine

	# To check if input files are there
	filename = 'ELT_2048_37m_11m_5mas_nospiders_cut.fits'
	my_file = str(conf['INPUT_DIR'] + filename)
	if (os.path.isfile(my_file) == False): 
		get_input_files()


def get_input_files():
	id2 = '1YR4G_8E7TpznTxumQA1zD4z_v4V5ZlF2'
	dir1 = '.'
	name  = 'fits_files.zip' 
	download_from_gdrive(id2, dir1 , name)
    
	with zipfile.ZipFile(name,'r') as zip_ref:
		zip_ref.extractall(conf['INPUT_DIR'])
	os.remove(name)





        





