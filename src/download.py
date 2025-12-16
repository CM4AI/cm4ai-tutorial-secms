import requests
from pyDataverse.api import NativeApi, DataAccessApi
import zipfile
import os
from pathlib import Path


cm4ai_doi = "doi:10.18130/V3/K7TGEM"
base_url = "https://dataverse.lib.virginia.edu"
target_filenames = [
    "cm4ai_mass-spec_KOLF2.zip",
    "cm4ai_mass-spec_MDA-MB-468.zip"
]

output_dir = "data/raw"

os.makedirs(output_dir, exist_ok=True)

api = NativeApi(base_url)

print("Downloading data set (this may take a while)...")
cm4ai_dataset = api.get_dataset(cm4ai_doi)
files = cm4ai_dataset.json()['data']['latestVersion']['files']

file_id = None
for file in files:
    if file['dataFile']['filename'] not in target_filenames:
        continue

    file_id = file['dataFile']['id']
    filename = file['dataFile']['filename']
    output_path = os.path.join(output_dir,filename)

    if not os.path.exists(output_path):
        download_url = f"{base_url}/api/access/datafile/{file_id}"
        response = requests.get(download_url)
        if response.status_code == 200:
            with open(output_path, 'wb') as f:
                f.write(response.content)
            print(f"Downloaded: {filename}")
        else:
            print(f"Failed to download file. Status code: {response.status_code}")
else:
    print("Target file not found in dataset.")

print("Download complete, extracing archives...")
for filename in os.listdir(output_dir):
    if filename.endswith('.zip'):
        zip_path = os.path.join(output_dir, filename)
        extract_path = os.path.join(output_dir, Path(zip_path).stem)
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(extract_path)
        os.remove(zip_path)

print(f"Complete! SEC-MS data extracted to: {output_dir}")
