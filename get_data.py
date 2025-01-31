import os
import requests
from subprocess import call

# Define the input file and output directory
input_file = '/home/shuaiw/Methy/benchmark.data.source.txt'
output_dir = '/home/shuaiw/methylation/data/published_data/fanggang/bax'
bam_dir = '/home/shuaiw/methylation/data/published_data/fanggang/bam'

# Create the output directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Function to download a file from a URL
def download_file(url, output_path):
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=1024):
                f.write(chunk)
        print(f"Downloaded {output_path}")
    else:
        print(f"Failed to download {url}")

# Read the input file and process each URL
annotation = None
i = 0
with open(input_file, 'r') as f:
    for line in f:
        if line.startswith('###'):
            annotation = line.strip().split('###')[1].strip()
            strain = annotation.split(',')[0].strip().replace(' ', '_').replace(';', '_').replace('(', '_').replace(')', '_')
            print (strain)
        elif line.startswith('http'):
            url = line.strip()
            file_name = os.path.basename(url)
            output_path = os.path.join(output_dir, file_name)
            url2 = url.replace('_p0.bas', '_p0.1.bax')
            file_name2 = os.path.basename(url2)
            output_path2 = os.path.join(output_dir, file_name2)
            url3 = url.replace('_p0.bas', '_p0.2.bax')
            file_name3 = os.path.basename(url3)
            output_path3 = os.path.join(output_dir, file_name3)
            url4 = url.replace('_p0.bas', '_p0.3.bax')
            file_name4 = os.path.basename(url4)
            output_path4 = os.path.join(output_dir, file_name4)


            download_file(url, output_path)
            download_file(url2, output_path2)
            download_file(url3, output_path3)
            download_file(url4, output_path4)

            cmd = f"""
            bax2bam -o {bam_dir}/{strain} {output_dir}/{file_name2} {output_dir}/{file_name3} {output_dir}/{file_name4}
            """ 
            print(cmd)
            os.system(cmd)
        i += 1
        # if i > 3:
        #     break

print("All files downloaded and converted.")