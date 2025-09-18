import tarfile
import io
from astropy.io import fits


def read_compressed_fits(compressed_file_path):
    # Open the tar.gz file
    """
    Read a FITS file compressed in a tar.gz archive and return the data as a numpy array.

    Parameters
    ----------
    compressed_file_path : str
        The path to the tar.gz file containing the FITS file.

    Returns
    -------
    data : numpy array
        The data in the FITS file as a numpy array.
    """
    with tarfile.open(compressed_file_path, 'r:gz') as tar:
        # Iterate over the members of the tar file
        for member in tar.getmembers():
            # Check if the member is a file (not a directory)
            if member.isfile():
                # Extract the file to a file-like object in memory
                file_obj = tar.extractfile(member)
                # Use astropy.io.fits to open the file-like object
                data = fits.getdata(io.BytesIO(file_obj.read()))
                # with fits.open(io.BytesIO(file_obj.read())) as hdul:
                #     # Access the data in the FITS file
                #     data = hdul[0].data  # Assuming the data is in the primary HDU
                #     header = hdul[0].header
                #     # Print or process the data and header as needed
                #     print(header)
                #     print(data)
    return data