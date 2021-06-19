"""
AUTHOR: Paarangat Pushkarna
DATE: 26/05/2021
MODIFIED: 30/05/2021
"""
from fitting import auto_fit
import os

RES = int(1e4)
WINDOW = 2

def main():
    work_path = os.getcwd()
    to_reconstruct = os.path.join(work_path, "to_reconstruct_kbeta.csv")
    ###USER_CHOICE###
    ct_folder_name = "sc_kb_ct_5s" ###Choose directory storing .ct files
    ######
    cts = os.path.join(work_path, ct_folder_name)

    auto_fit(cts, to_reconstruct, RES, ct_folder_name, WINDOW)

if __name__ == "__main__":
    main()
