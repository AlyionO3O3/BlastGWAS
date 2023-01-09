import pandas as pd
import requests
import wget, tarfile
import os
import gzip
import coreapi
import numpy as np
from liftover import ChainFile

class Rmcga():
    def __init__(self):
        self.DEG = None
        self.bim = None
        self.refgene = None
        self.tf_list = None
        self.tfYesmotif_dict = {}
        self.tfNomotif_list = []
        self.load_deg()
        self.load_refgene()
        self.load_gwas()
        self.load_tf()
        self.scrape_jaspar()
        self.paste_chip_url()
        self.lift_beds()

#load DEG file
    def load_deg(self):
        try:
            print("Loading DEG list file...")
            self.DEG = pd.read_csv("DEG_list.csv", sep=",")
            print("DEG list file loaded...")
        except:
            print("Fail to load DEG list file...")

    def load_gwas(self):#load GWAS pre QC bim file to acquire SNP position
        try:
            print("Loading GWAS bim file...")
            self.bim = pd.read_csv("combined.TWB1.TWB2.high.confidence.v3.bim", sep="\t")
            print("GWAS bim file loaded...")
        except:
            print("Fail to load GWAS bim file...")

    def load_refgene(self,url='https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz'):#load UCSC gene anotation file
        try:
            print("Pulling refgene request...")
            request = requests.get(url, allow_redirects=True)
            file_name = url.split('/')[-1]
            print("Writing refgene file...")
            open(file_name, "wb").write(request.content)
            self.refgene = pd.read_csv(file_name, compression='gzip', header=None, sep='\t', quotechar='"',
                        on_bad_lines='skip')
            os.remove(file_name)
            print("Refgene file created...")
        except:
            print("Fail to create refgene file...")

    def load_tf(self):
        try:
            print("Loading TF list...")
            self.tf_list = []
            with open("254_TFs.txt") as f:
                for line in f.readlines():
                    self.tf_list.append(line.rsplit())
            self.tf_list = [x for item in self.tf_list for x in item]
            print("TF list loaded...")
        except:
            print("Fail to load TF list...")

    def get_matrix_id(self,tf):
        # Initialize a client & load the schema document
        client = coreapi.Client()
        schema = client.get("https://jaspar.genereg.net/api/v1/docs/")
        # Interact with the API endpoint
        action = ["matrix", "list"]
        params = {
            "search": tf,
            "tax_id": "9606",
            "version": 'latest',
            "data_type": "ChIP-seq",
            "release": "2022"
        }
        result = client.action(schema, action, params=params)
        if tf.lower() in result['results'][0]['name'].lower() or result['results'][0]['name'].lower() in tf.lower():
            return result['results'][0]['matrix_id']

    def scrape_jaspar(self):
        print("Start JASPAR scraping...")
        i = len(self.tf_list)
        for tfs in self.tf_list:
            try:
                MOTIF = get_matrix_id(tfs)
                if MOTIF == None:
                    print(f'Can not find JASPAR MOTIF ID for {tfs}!')
                else:
                    self.tfYesmotif_dict.setdefault(tfs, MOTIF)
                    print(f'The JASPAR MOTIF ID for {tfs} is {MOTIF}!')
            except:
                self.tfNomotif_list.append(tfs)
                print(f'Can not find JASPAR MOTIF ID for {tfs}!')
            i -= 1
            if i != 0:
                print(f'{i} TFs left!')
            else:
                print('Scraping done!!')
            print(f'We found {len(self.tfYesmotif_dict)} MOTIFs from our uploaded TFs by scraping JASPAR!')
            print(f'{len(self.tfNomotif_list)} other TFs failed, listed below:')
            print(self.tfNomotif_list)
    def paste_chip_url(self):
        print("Pasting JASPAR url...")
        self.tfmotifURL_dict = {}
        for i in self.tfYesmotif_dict.keys():
            motif_id = self.tfYesmotif_dict[i]
            motif_url = 'https://jaspar.genereg.net/download/data/2022/bed/' + motif_id + '.bed'
            self.tfmotifURL_dict.setdefault(i, motif_url)
        print("Done Pasting JASPAR url...")

    def load_chip_file(self,url):
        request = requests.get(url, allow_redirects=True)
        file_name = url.split('/')[-1]
        open(file_name, "wb").write(request.content)
        ChIP = pd.read_csv(file_name, header=None, sep='\t', quotechar='"', on_bad_lines='skip')
        os.remove(file_name)
        return ChIP
    def liftover(self,chrom,pos):
        converter = ChainFile('hg19ToHg38.over.chain.gz', 'hg18', 'hg38')
        return converter[chrom][pos]
    def lift_beds(self):
        for i in self.tfmotifURL_dict.value:
            ChIP = self.load_chip_file(i)
            ChIP.columns = ["CHR","Start","End","POS","Peak","Strand"]
            if ChIP.POS


    def mapping_tfbs(self):
        pass

    def gene_annotion(self):
        pass

    def dgs_annotion(self):
        pass

if __name__ == "__main__": main()