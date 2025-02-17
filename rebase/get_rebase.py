from Bio import Entrez
import xmltodict
import re
from collections import defaultdict
## 
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager
from selenium.webdriver.common.by import By
import time
import pandas as pd


# Set NCBI Entrez email (required)
Entrez.email = "wshuai294@gmail.com"  # Replace with your email

def get_sra_ids_from_biosample(biosample_id):
    """Searches NCBI for a BioSample ID and retrieves all linked SRA Run IDs."""
    try:
        # Step 1: Search for the BioSample ID in the NCBI database
        search_handle = Entrez.esearch(db="biosample", term=biosample_id)
        search_results = Entrez.read(search_handle)
        search_handle.close()

        if not search_results["IdList"]:
            print(f"No BioSample entry found for ID: {biosample_id}")
            return None

        biosample_ncbi_id = search_results["IdList"][0]  # Get the first matching BioSample entry

        # Step 2: Fetch the full BioSample metadata
        fetch_handle = Entrez.efetch(db="biosample", id=biosample_ncbi_id, retmode="xml")
        biosample_data = xmltodict.parse(fetch_handle.read())
        fetch_handle.close()

        # Step 3: Extract all SRA Run IDs from BioSample metadata
        sra_ids = []
        links = biosample_data["BioSampleSet"]["BioSample"].get("Links", {}).get("Link", [])

        if isinstance(links, list):  # If multiple links exist
            for link in links:
                if "sra" in link["@type"].lower():
                    sra_ids.append(link["@target"])
        elif isinstance(links, dict):  # If there's only one link
            if "sra" in links["@type"].lower():
                sra_ids.append(links["@target"])

        if not sra_ids:
            print(f"No direct SRA links found for BioSample {biosample_id}, checking SRA database...")

            # Step 4: If no direct SRA links, search in the SRA database
            sra_handle = Entrez.esearch(db="sra", term=biosample_id)
            sra_results = Entrez.read(sra_handle)
            sra_handle.close()

            if sra_results["IdList"]:
                sra_ids = sra_results["IdList"]

        if sra_ids:
            print(f"BioSample {biosample_id} -> SRA Runs: {sra_ids}")
            return sra_ids
        else:
            print(f"No SRA Runs found for BioSample {biosample_id}.")
            return None

    except Exception as e:
        print(f"Error retrieving SRA Runs: {e}")
        return None



def convert_sra_uid_to_srr(sra_uid):
    """Convert an SRA UID to its standard SRA Run ID (SRRxxxxxxx)."""
    try:
        # Fetch metadata for the given SRA UID
        handle = Entrez.efetch(db="sra", id=sra_uid, retmode="xml")
        sra_metadata = xmltodict.parse(handle.read())
        handle.close()

        # Extract the SRA Run ID (SRRxxxxxxx)
        run_id = sra_metadata["EXPERIMENT_PACKAGE_SET"]["EXPERIMENT_PACKAGE"]["RUN_SET"]["RUN"]["@accession"]

        print(f"SRA UID {sra_uid} -> Standard SRA ID: {run_id}")
        return run_id

    except Exception as e:
        print(f"Error converting SRA UID {sra_uid}: {e}")
        return None
    


def get_http_links_from_sra_page(sra_id):
    """Retrieve all HTTP links from NCBI SRA Run Browser for a given SRA Run ID."""
    try:
        # Construct the URL for the SRA Run Browser
        url = f"https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc={sra_id}&display=data-access"
        print(f"Opening page: {url}")

        # Set up Selenium WebDriver
        options = webdriver.ChromeOptions()
        options.add_argument("--headless")  # Run Chrome in headless mode
        options.add_argument("--no-sandbox")
        options.add_argument("--disable-dev-shm-usage")

        service = Service(ChromeDriverManager().install())
        driver = webdriver.Chrome(service=service, options=options)

        # Open the webpage
        driver.get(url)
        time.sleep(5)  # Allow JavaScript to load the page

        # Extract all HTTP links
        http_links = set()
        elements = driver.find_elements(By.TAG_NAME, "a")

        for element in elements:
            href = element.get_attribute("href")
            if href and href.startswith("http"):
                http_links.add(href)

        driver.quit()

        if http_links:
            print(f"Found {len(http_links)} links for {sra_id}:")
            bam_link = []
            for link in sorted(http_links):
                if re.search("sra-pub-src", link):
                    print(link)
                    bam_link.append(link)
            return bam_link
        else:
            print(f"No HTTP links found for {sra_id}.")
            return None

    except Exception as e:
        print(f"Error retrieving HTTP links for {sra_id}: {e}")
        return None

def each_biosample(biosample_id):
    # biosample_id = "SAMN04419113"
    uid_list = get_sra_ids_from_biosample(biosample_id)

    ## if uid_list is None
    if uid_list is None:
        return {}


    biosample_dict = defaultdict(dict)
    ## transfer SRA UID to SRR
    standard_sra_ids = []
    for uid in uid_list:
        sra_id = convert_sra_uid_to_srr(uid)

        # sra_id = "SRR5937949"  # Replace with actual SRA Run ID
        bam_link = get_http_links_from_sra_page(sra_id)
        biosample_dict[sra_id] = bam_link
        standard_sra_ids.append(sra_id)

    # print(standard_sra_ids)
    # print(biosample_dict)
    return biosample_dict


if __name__ == "__main__":
    # Example usage
    # biosample_id = "SAMN07447437"  # Replace with your BioSample ID
    df = pd.read_csv("basemodification.csv")
    biosample_id_uniq = set()
    data = []
    for index, row in df.iterrows():
        biosample_id = row["BioSample"]
        
        if biosample_id in biosample_id_uniq:
            continue
        print (biosample_id)
        biosample_id_uniq.add(biosample_id)

        biosample_dict = each_biosample(biosample_id)
        
        for sra in biosample_dict:
            if biosample_dict[sra] is None:
                continue
            for bam in biosample_dict[sra]:
                data.append([row['BioProject'],row['Organism'], biosample_id, sra, bam])
                print (row['BioProject'],row['Organism'], biosample_id, sra, bam)
        print("Done", len(data))
        # if len(data) > 10:
        #     break
    df = pd.DataFrame(data, columns=["BioProject", "Organism", "BioSample", "SRA", "BAM"])
    df.to_csv("rebase_bam_data.csv", index=False)



