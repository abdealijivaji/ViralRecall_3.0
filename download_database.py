from pathlib import Path
import requests



def main() :
    base_dir = Path(__file__).parent
    # hmm_dir = (base_dir / 'hmm').mkdir(exist_ok=)
    dbsrc = "https://zenodo.org/records/12666277/files/hmm.tar.gz"

    requests.get(dbsrc)



    pass




if __name__ == "__main__" :
    main()
