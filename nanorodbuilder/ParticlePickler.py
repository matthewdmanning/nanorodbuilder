import pickle
import NanoParticle
import pathlib

pickle_dir = '/home/mdmannin/mnt/storage/nanoparticles/python_objects'

def pickle_nanoparticle(nanoparticle_object, nanoparticle_name):
    np_file_name = '{}/{}.pickle'.format(pickle_dir, nanoparticle_name)
    if not pathlib.Path(np_file_name).is_file():
        np_file = open(np_file_name, 'wb')
        pickle.dump(nanoparticle_object, np_file)
        print('File pickled: {}'.format(np_file_name))
        np_file.close()
        return True
    else:
        print('File exists, not pickling: {}'.format(np_file_name))
        return False

def unpickle_nanoparticle(nanoparticle_name):
    np_file_name = '{}/{}.pickle'.format(pickle_dir, nanoparticle_name)
    if pathlib.Path(np_file_name).is_file():
        np_file = open(np_file_name, 'r')
        nanoparticle = pickle.load(np_file)
        print('Nanoparticle loaded from pickled file: {}'.format(np_file_name))
        np_file.close()
        return nanoparticle
    else:
        print('Nanoparticle pickle file not found.')
        return None