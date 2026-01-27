import Bio
from Bio import SeqIO
from pathlib import Path
import pandas as pd

gbkFile = Path("/home/abdeali/viralR_test_output/Chlamy_punui/Chlamy_punui_contig.gbk")

gbk = SeqIO.read(open(gbkFile, "r"), "genbank")

# print(f"IDs in file: {gbk.features[2]} ")

def features_to_dataframe(recs, cds=False):
    """Get genome records from a biopython features object into a dataframe
      returns a dataframe with a row for each cds/entry"""

    genome = recs[0]
    #preprocess features
    allfeat = []
    for (item, f) in enumerate(genome.features):
        x = f.__dict__
        q = f.qualifiers
        x.update(q)
        d = {}
        d['start'] = f.location.start
        d['end'] = f.location.end
        d['strand'] = f.location.strand
        for i in f.keys:
            if i in x:
                if type(x[i]) is list:
                    d[i] = x[i][0]
                else:
                    d[i] = x[i]
        allfeat.append(d)

    import pandas as pd
    df = pd.DataFrame(allfeat,columns=f.keys)
    df['length'] = df.translation.astype('str').str.len()
    #print (df)
    # df = check_tags(df)
    # if cds == True:
    #     df = get_cds(df)
    #     df['order'] = range(1,len(df)+1)
    #print (df)
    if len(df) == 0:
        print ('ERROR: genbank file return empty data, check that the file contains protein sequences '\
               'in the translation qualifier of each protein feature.' )
    return df

features = features_to_dataframe(gbk)

print(features)