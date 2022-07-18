import funsies as f
import argparse, os, json

parser = argparse.ArgumentParser(description="Optimize Lewis acid Epoxide complex")

parser.add_argument("-SMILES", "--smiles", help="smiles", type=str)
parser.add_argument("-c", "--charge", help="charge", type=int, default=0)
parser.add_argument("-s", "--spin", help="spin", type=int, default=1)
parser.add_argument("-n", "--ntasks", help="number of cores", type=str, default="4")

args = parser.parse_args()

SMILES = args.smiles

complex_ID_dict = json.load(open("LAs/complex_id.json"))
ID = "_".join(complex_ID_dict[SMILES].split())

obabel_geom = ID + "_" + "obabel.xyz"
xtb_opt_out = ID + "_" + "xtb_opt.xyz"

with f.Fun():
    
    obabel = f.shell(
        f"obabel -:'{SMILES}' --gen3d -O obabel.xyz",
        out=["obabel.xyz"],
        )

    xtb_opt = f.shell(
        f"xtb input.xyz --opt --alpb acetonitrile -P 4",
        inp={"input.xyz": obabel.out["obabel.xyz"]},
        out=["xtbopt.xyz"],
        env={"OMP_NUM_THREADS": str(args.ntasks)},
        )

    f.execute(xtb_opt)
    f.wait_for(xtb_opt)

    os.makedirs("xtb_opts", exist_ok=True)
    f.takeout(xtb_opt.out["xtbopt.xyz"], "xtb_opts/"+xtb_opt_out)
    




