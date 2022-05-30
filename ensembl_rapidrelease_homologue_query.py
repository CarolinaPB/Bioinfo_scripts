# selenium 4
import time
from selenium import webdriver
from selenium.webdriver.firefox.service import Service
from webdriver_manager.firefox import GeckoDriverManager
from bs4 import BeautifulSoup
from selenium.webdriver import FirefoxOptions
import pandas as pd
from ratelimit import limits, sleep_and_retry
from getpass import getpass
import os
import argparse

### USAGE ###
# ensembl_rapidrelease_homologue_query.py -i <ids.txt> -u <https://rapid.ensembl.org/species_code/Gene/Compara_Homolog?g=>
# -i is a file with Ensembl gene IDs for your species' genes that you want to find homologues for
# -u is teh rapid release homologue page url, something like this https://rapid.ensembl.org/Meleagris_gallopavo_GCA_905368555.1/Gene/Compara_Homolog?g=

# Check below for other fields you should change

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--ids", help = "Text file containing one Ensembl gene id per line")
parser.add_argument("-u", "--url", help = "Rapid release homologue URL")

args = parser.parse_args()

os.environ['GH_TOKEN'] = getpass("Github token: ")

one_hour = 3600

@sleep_and_retry
@limits(calls=5000 , period=one_hour)
def parse(url):
    opts = FirefoxOptions()
    opts.add_argument("--headless")
    response = webdriver.Firefox(service=Service(GeckoDriverManager().install()),options=opts)
    response.get(url)
    
    print(url)
    time.sleep(5)
    sourceCode=response.page_source
    response.quit() 
    return  sourceCode

with open(args.ids) as f:
    gene_list = f.read().splitlines()

# TODO: change species names, add or remove items from the list
final_df = pd.DataFrame(index = gene_list, columns=["Mgal_5.1_ID", "Chicken_ID"])

if not os.path.isfile("processed_genes.txt"):
    with open("processed_genes.txt", "w") as gene_log:
        pass

for gene, row in final_df.iterrows():
    with open('processed_genes.txt') as f:
        if gene in f.read():
            print(f"Gene {gene} already processed")
        else:
            url = f"{args.url}{gene}"

            soup = BeautifulSoup(parse(url), "lxml")

            if not soup.body.findAll(text='No homologues have been found for this gene.'):
                x = soup.find_all("table")

                my_table = x[1]

                table_body = my_table.find("tbody")

                data = []

                rows = table_body.find_all('tr')
                for row in rows:
                    cols = row.find_all('td')
                    cols = [ele.text.strip() for ele in cols]
                    data.append([ele for ele in cols if ele]) # Get rid of empty values

# TODO: change species names, add or remove items if statements
# Turkey and Chicken correspond to species in the homologue gene table, change to your species' name
                for el in data:
                    if "Turkey" in el:                
                        final_df.loc[gene, "Mgal_5.1_ID"] = el[1] if el[1].startswith("ENS") else el[2]

                    elif "Chicken" in el:
                        final_df.loc[gene, "Chicken_ID"] = el[1] if el[1].startswith("ENS") else el[2]
            else:
                print(f"No homologue for gene: {gene}")
                final_df.loc[gene, "Mgal_5.1_ID"] = None
                final_df.loc[gene, "Chicken_ID"] = None

            with open("processed_genes.txt", "a") as f:
                f.write(gene + "\n")

            final_df.to_csv("homologue_table.csv", mode='a')






