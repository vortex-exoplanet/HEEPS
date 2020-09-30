import pickle
import os.path

def save2pkl(data, name, dir_output='output_files', band='L', mode='RAVC', **conf):

    ''' Save data to a pickle file. '''

    filename = os.path.join(dir_output, '%s_%s_%s.pkl'%(name, band, mode))
    with open(filename, 'wb') as file:
        pickle.dump(conf, file)
        file.close()

    return filename