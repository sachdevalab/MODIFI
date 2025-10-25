import requests
from bs4 import BeautifulSoup

def get_genome_links():
    """Scrape the NEB Genomes page and extract genome download links into a dictionary."""
    base_url = "https://tools.neb.com/genomes/"
    genome_dict = {}

    try:
        response = requests.get(base_url)
        if response.status_code != 200:
            print(f"Failed to access {base_url}")
            return None
        
        soup = BeautifulSoup(response.text, "html.parser")

        # Find all links containing genome data
        for link in soup.find_all("a", href=True):
            href = link["href"]
            # if href.endswith(".fasta") or href.endswith(".gz"):  # Adjust for different genome file formats
            full_url = base_url + href if href.startswith("/") else href
            print (full_url)
            genome_name = href.split("/")[-1]  # Extract filename as genome name
            genome_dict[genome_name] = full_url

        if genome_dict:
            print(f"Found {len(genome_dict)} genome files:")
            for genome, url in genome_dict.items():
                print(f"{genome}: {url}")
        else:
            print("No genome files found.")

        return genome_dict

    except Exception as e:
        print(f"Error scraping genome links: {e}")
        return None

# Run the script
genome_links_dict = get_genome_links()

