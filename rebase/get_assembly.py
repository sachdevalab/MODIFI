from Bio import Entrez
import xmltodict
import requests
import re

# Set your NCBI Entrez email (required)
Entrez.email = "your_email@example.com"  # Replace with your email

def get_assembly_from_biosample(biosample_id):
    """Retrieve the Assembly download link and assembly level for a given BioSample ID."""
    try:
        # Step 1: Search for the Assembly linked to the BioSample
        handle = Entrez.esearch(db="assembly", term=biosample_id)
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            print(f"No Assembly found for BioSample ID: {biosample_id}")
            return None

        assembly_id = record["IdList"][0]

        # Step 2: Fetch the Assembly metadata
        fetch_handle = Entrez.efetch(db="assembly", id=assembly_id, retmode="xml")
        assembly_data = xmltodict.parse(fetch_handle.read())
        fetch_handle.close()

        # Step 3: Extract Assembly details
        assembly_entry = assembly_data["AssemblySet"]["Assembly"]
        ftp_path = assembly_entry["FtpPath_RefSeq"]  # Use RefSeq FTP link (preferred)
        if not ftp_path:
            ftp_path = assembly_entry["FtpPath_GenBank"]  # Fallback to GenBank FTP link

        assembly_level = assembly_entry["AssemblyStatus"]  # Assembly level (Complete, Scaffold, Contig, etc.)

        # Construct full FASTA download link
        fasta_link = f"{ftp_path}/{ftp_path.split('/')[-1]}_genomic.fna.gz"

        print(f"BioSample {biosample_id} -> Assembly Level: {assembly_level}")
        print(f"Download Link: {fasta_link}")

        return {"assembly_level": assembly_level, "download_link": fasta_link}

    except Exception as e:
        print(f"Error retrieving Assembly for BioSample {biosample_id}: {e}")
        return None

# Example usage
biosample_id = "SAMN12345678"  # Replace with your BioSample ID
get_assembly_from_biosample(biosample_id)

