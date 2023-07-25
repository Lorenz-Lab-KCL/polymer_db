import os
import json
import shutil

from polymer_db.orca_pella_db import OrcaPellaReader


if __name__ == "__main__":
    target_dir = "poly_json_large_37292_incomplete_smiles"
    os.mkdir(target_dir)
    os.mkdir("xyz")

    try:
        filenames_schemas = OrcaPellaReader.get_all_jsons_from_folder(
            folder_path="poly_db_large_37292",
            with_multiprocessing=True,
            track_progress=True,
        )
        for filename_schema in filenames_schemas:
            if filename_schema is not None:
                filename, schema = filename_schema
                file_id = filename.split("/")[1]
                with open(f"{target_dir}/{file_id}.json", "w") as f:
                    json.dump(schema, f, indent=4)
    finally:
        shutil.rmtree("xyz")
