
from pathlib import Path
import requests, time, progressbar, shutil
import numpy as np
from src.utils import prep_hmm
import hashlib, argparse, os

def check_file_hash(filepath : Path, expected_hash : str) -> bool :
    ''' Check the md5 hash of a file to ensure it matches the expected hash '''
    
    with open(filepath, "rb") as f:
        file_hash = hashlib.md5(f.read()).hexdigest()
    return file_hash == expected_hash

def download_file(url : str, filepath : Path) -> None:
    compressed_hash = "fa81765a316c5d84c67c9da12a10ee9e"
    if filepath.exists() and check_file_hash(filepath, compressed_hash) :
        print(f"Database is already downloaded")
        return
    n_chunk=1000
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
            time.sleep(0.05)
    if filepath.exists() and check_file_hash(filepath, compressed_hash) :
        print(f"Finished downloading HMM database")
    else :
        print(f"Error in downloading database. Please try again.")
    return

def parse_args(argv=None) :
    args_parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, 
                    description="Download database for Viralrecall v3.0", 
                    epilog='*******************************************************************\n\n*******************************************************************')
    args_parser.add_argument('-d', '--db_dir', required=False, default=Path.cwd(),
                             help='Directory to download the HMM database. Default is current working directory')
    args_parser.add_argument('-n', '--name' , required=False, default='hmm',
                             help='Name to give the database directory. Default is "hmm"')
    args_parser = args_parser.parse_args()
    
    return args_parser


    # return outname.with_suffix('.h3m')

def main() :
    args_list = parse_args()
    out_dir = Path(args_list.db_dir).expanduser()
    project_name = args_list.name 
    
    print(f"Downloading database to {out_dir}")
    dbsrc = "https://zenodo.org/records/17859729/files/hmm.tar.gz"
    
    db_dir = out_dir / project_name
    db_file = db_dir.with_suffix('.tar.gz') 
    download_file(dbsrc, db_file)

    if not db_dir.is_dir() :
        shutil.unpack_archive(db_file, out_dir)
        print("Finished unpacking")
    os.remove(db_file)
    print(f"Preparing HMM database: {db_dir}")
    prep_hmm(db_dir)



if __name__ == "__main__" :
    main()
