"""Figure 3 paired patient-level summary using the final 87 CL/CM pairs."""
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import wilcoxon
FEATURES=['Grad_Angle_DZ_SZ','Surface_Roughness','SZ_ChondronDensity_mm2','Median_Thickness_um']
SCORE_COLS=['True_Total_OARSI_Mean','True_Total_HHGS_Mean','True_Structure_Mean']
def paired_summary(df):
    if df.duplicated(['Patient','Region']).any():raise ValueError('Duplicate patient-region records')
    records=[];pairs=[]
    for feature in FEATURES+SCORE_COLS:
        missing=set(['Patient','Region',feature])-set(df.columns)
        if missing:raise KeyError(missing)
        tab=df.pivot(index='Patient',columns='Region',values=feature)
        pair=tab[['Lateral','Medial']].dropna()
        if len(pair) != 87:
            raise ValueError(
                f"{feature}: expected 87 complete CL/CM pairs; found {len(pair)}"
            )
        # Original paired Wilcoxon method and two-sided null; no statistical value is hard-coded.
        stat,pval=wilcoxon(pair['Lateral'].to_numpy(float),pair['Medial'].to_numpy(float),alternative='two-sided',zero_method='wilcox')
        records.append(dict(Feature=feature,N_pairs=len(pair),Lateral_mean=pair['Lateral'].mean(),Medial_mean=pair['Medial'].mean(),Wilcoxon_stat=stat,Wilcoxon_P=pval))
        q=pair.reset_index().rename(columns={'Lateral':'CL','Medial':'CM'});q.insert(0,'Feature',feature);pairs.append(q)
    return pd.DataFrame(records),pd.concat(pairs,ignore_index=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--merged-csv',required=True);ap.add_argument('--output-dir',default='outputs/simtk');args=ap.parse_args()
    d=pd.read_csv(args.merged_csv);o=Path(args.output_dir);o.mkdir(parents=True,exist_ok=True)
    summ,pairs=paired_summary(d);summ.to_csv(o/'Figure3_Paired_Summary.csv',index=False);pairs.to_csv(o/'Figure3_Paired_SourceData.csv',index=False)
    print(summ.to_string(index=False))
