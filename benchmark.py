from utils import *
import pandas as pd
dft_dict = {'HF': ade.OptKeywords(['HF', 'DEF2-TZVPP', 'DEF2-TZVPP/C', 'def2/J','RIJCOSX','CPCM']),
                'RI-MP2': ade.OptKeywords(['RI-MP2', 'DEF2-TZVPP', 'DEF2-TZVPP/C', 'def2/J','RIJCOSX','CPCM']),
                'BP': ade.OptKeywords(['BP', 'DEF2-TZVPP', 'DEF2-TZVPP/C', 'def2/J','RIJCOSX','CPCM']),
                'PBEh-3c': ade.OptKeywords(['PBEh-3c', 'DEF2-TZVPP', 'DEF2-TZVPP/C', 'def2/J','RIJCOSX','CPCM']),
                'B97': ade.OptKeywords(['B97', 'DEF2-TZVPP', 'DEF2-TZVPP/C', 'def2/J','RIJCOSX','CPCM']),
                'BLYP': ade.OptKeywords(['BLYP', 'DEF2-TZVPP', 'DEF2-TZVPP/C', 'def2/J','RIJCOSX','CPCM']),
                'PBE': ade.OptKeywords(['PBE', 'DEF2-TZVPP', 'DEF2-TZVPP/C', 'def2/J','RIJCOSX','CPCM', 'D4']),
                'PBE0': ade.OptKeywords(['PBE0', 'DEF2-TZVPP', 'DEF2-TZVPP/C', 'def2/J','RIJCOSX','CPCM', 'D4']),
                'B3LYP': ade.OptKeywords(['B3LYP', 'DEF2-TZVPP', 'DEF2-TZVPP/C', 'def2/J','RIJCOSX','CPCM', 'D4']),
                'M062X': ade.OptKeywords(['M062X', 'DEF2-TZVPP', 'DEF2-TZVPP/C', 'def2/J','RIJCOSX','CPCM', 'D4']),
                'M06L': ade.OptKeywords(['M06L', 'DEF2-TZVPP', 'DEF2-TZVPP/C', 'def2/J','RIJCOSX','CPCM', 'D4']),
                'TPSSh': ade.OptKeywords(['TPSSh', 'DEF2-TZVPP', 'DEF2-TZVPP/C', 'def2/J','RIJCOSX','CPCM']),
                'wB97': ade.OptKeywords(['wB97', 'DEF2-TZVPP', 'DEF2-TZVPP/C', 'def2/J','RIJCOSX','CPCM', 'D4']),
                'wB97X': ade.OptKeywords(['wB97X', 'DEF2-TZVPP', 'DEF2-TZVPP/orcaC', 'def2/J','RIJCOSX','CPCM', 'D4']),
                'wB97X-V': ade.OptKeywords(['wB97X-V', 'DEF2-TZVPP', 'DEF2-TZVPP/C', 'def2/J','RIJCOSX','CPCM']),
                'RI-B2PLYP': ade.OptKeywords(['RI-B2PLYP', 'DEF2-TZVPP', 'DEF2-TZVPP/C', 'def2/J','RIJCOSX','CPCM']),
                'DLPNO-CCSD(T)': ade.OptKeywords(['DLPNO-CCSD(T)', 'DEF2-TZVPP', 'DEF2-TZVPP/C', 'def2/J','RIJCOSX', 'CPCM']),
            }

df = pd.DataFrame()
methods=['PBE', 'PBE0', 'M062X', 'M06L', 'B3LYP', 'wB97', 'wB97X', 'DLPNO-CCSD(T)']
df["method"] = methods

df = benchmark_DFT(ID="boron_trifluoride", dft_dict=dft_dict, methods=methods, solvent="Acetonitrile", df=df)
df = benchmark_DFT(ID="aluminium_chloride", dft_dict=dft_dict, methods=methods, solvent="Acetonitrile", df=df)
print(df)

