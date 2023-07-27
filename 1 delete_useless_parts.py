import logging
import os
import shutil
from tqdm import tqdm


logging.basicConfig(format="%(levelname)s: %(message)s", filename="error raw.log")


def _get_num_at_end(_str: str, key: str = "_"):
    return _str.split(key)[-1]


def save_orca_files_with_id(orca_folder_path: str, target_root_dir: str):
    # get uuid
    orca_folder_path += "/orca"
    id_parts = orca_folder_path.split("/")
    file_id = []
    for id_part in id_parts[-4:-1]:
        file_id.append(_get_num_at_end(id_part))
    file_id = "_".join(file_id)

    # copy correct files to target folder
    chosen_files = []
    for file in os.listdir(orca_folder_path):
        if "log" in file or ("mol" in file and "smi" in file):
            chosen_files.append(file)

    if not len(chosen_files) == 2:
        logging.error(f"{file_id} folder is unnatural: {orca_folder_path}")
        return

    for chosen_file in chosen_files:
        source_file_path = os.path.join(orca_folder_path, chosen_file)
        target_file_path = os.path.join(target_root_dir, file_id, chosen_file)
        os.makedirs(os.path.dirname(target_file_path), exist_ok=True)
        shutil.copyfile(source_file_path, target_file_path)


if __name__ == "__main__":
    source_dir = "/Volumes/Seagate Portable Drive/polymer_database/res_ext_bridges_3"
    target_dir = "/Users/kouroshzarei/Desktop/poly_db"

    parts = [os.path.join(source_dir, part) for part in os.listdir(source_dir)]
    parts = [i for i in parts if os.path.isdir(i)]

    # print(bridges[:3])

    # parts = [
    #     os.path.join(bridge, part)
    #     for bridge in tqdm(bridges[:3])
    #     for part in os.listdir(bridge)
    # ]
    # parts = [i for i in parts if os.path.isdir(i)]

    sub_parts = [
        os.path.join(part, sub_part)
        for part in tqdm(parts)
        for sub_part in os.listdir(part)
    ]
    sub_parts = [
        i for i in tqdm(sub_parts) if os.path.isdir(i) if "orca" in os.listdir(i)
    ]

    for sub_part in tqdm(sub_parts):
        save_orca_files_with_id(sub_part, target_dir)
