
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages  
import pandas as pd
from pathlib import Path
import logging

plt.style.use('petroff10')
logging.getLogger('matplotlib.font_manager').disabled = True


def plot_vreg(annot_tbl : pd.DataFrame, minscore: int, 
              vreg_coords : dict, out_base: Path) -> None :
    out_file = out_base.parent / 'plots.pdf' 
    with PdfPages(out_file) as of :
        for key, coords in vreg_coords.items():
            df = annot_tbl.loc[annot_tbl['contig'] == key]
            x = df['pstart']/1_000_000
            y = df['rollscore']
            starts = [i / 1_000_000 for i in coords[0]]
            ends = [i / 1_000_000 for i in coords[1]]
            
            fig, ax = plt.subplots(layout="constrained", figsize=(10, 6))

            

            plt.plot(x, y, linewidth = '2')
            plt.axhline(minscore, color='red', linestyle='--', label='Min Score Threshold')
            plt.xlabel("Genome Position (Mb)")
            plt.ylabel("Score")
            plt.title(f"Viralrecall Score Plot for {key}")
            for i, _ in enumerate(starts):
                pos_annot = ((starts[i] + ends[i])/2 , max(y)*0.9)
                ax.annotate(f'vregion_{i+1}', xy = pos_annot, xytext=pos_annot, 
                            fontsize=10, ha='center')
                plt.fill_betweenx(y=[0, max(y)], x1=starts[i], x2=ends[i], 
                            color='green', alpha=0.1)
            of.savefig()
            plt.close()

