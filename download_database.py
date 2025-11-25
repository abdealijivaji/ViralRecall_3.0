from pathlib import Path
import requests, time, progressbar, shutil
import numpy as np
from src.utils import prep_hmm

def download_file(url : str, filepath : Path, n_chunk=1000):
    if filepath.exists() :
        print(f"Database is already downloaded")
        return
    r = requests.get(url, stream=True)
    
    # Estimates the number of bar updates
    block_size = 1024
    file_size = int(r.headers.get('Content-Length', 0))
    print(f"Downloading database of total size: {int(file_size/(block_size**2))} MB")
    num_bars = np.ceil(file_size / (n_chunk * block_size))
    bar =  progressbar.ProgressBar(maxval=num_bars).start()
    with open(filepath, 'wb') as f:
        for i, chunk in enumerate(r.iter_content(chunk_size=n_chunk * block_size)):
            f.write(chunk)
            bar.update(i+1)
            # Add a little sleep so you can see the bar progress
            time.sleep(0.05)
    print(f"\nFinished downloading HMM database")
    return



    # return outname.with_suffix('.h3m')

def main() :
    base_dir = Path("~").expanduser() # Path(__file__).parent
    db_dir = base_dir / 'hmm'
    dbsrc = "https://zenodo.org/records/12666277/files/hmm.tar.gz"
    db_file = db_dir.with_suffix('.tar.gz')
    download_file(dbsrc, db_file)
    
    if not db_dir.is_dir() :
        shutil.unpack_archive(db_file, base_dir)
        print("Finished unpacking")
    
    prep_hmm(db_dir)





if __name__ == "__main__" :
    main()
