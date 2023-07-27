import os
import tarfile
from tqdm import tqdm


def make_tarfile(source_dir, output_filename):
    with tarfile.open(output_filename, "w:gz") as tar:
        # Get the list of files in the source directory
        files = []
        for root, _, filenames in os.walk(source_dir):
            for filename in filenames:
                files.append(os.path.join(root, filename))

        # Add files to the tar archive with progress bar
        progress_bar = tqdm(total=len(files), desc="Creating tar archive", unit="file")
        for file in files:
            tar.add(file, arcname=os.path.relpath(file, source_dir))
            progress_bar.update(1)
        progress_bar.close()


if __name__ == "__main__":
    fol_path = "poly_json_small_1122_incomplete_smiles"
    tar_path = "poly_comp_json_small_1122_incomplete_smiles.tar.gz"
    make_tarfile(fol_path, tar_path)
