import numpy as np
import pandas as pd
import requests
import os
import coreapi
import subprocess
from liftover import ChainFile
from tqdm import tqdm
import multiprocessing as mp
import time
start_time = time.time()
tqdm.pandas()
# python3 -m venv tutorial-env
# source tutorial-env/bin/activate
class Rmcga():
    def __init__(self):
        self.DEG = None
        self.gwas = None
        self.refgene = None
        self.gwas_result = None
        self.tf_list = []
        self.tfyesmotif_dict = {}
        self.tfnomotif_list = []
        self.tfmotifURL_dict = None

    def totxt(self, df, name, path=os.getcwd() + "/", h=True):
        chunks = np.array_split(df.index, 100)
        for chunck, subset in enumerate(tqdm(chunks)):
            if chunck == 0:  # first row
                df.loc[subset].to_csv(path + name, mode='w', header=h, index=False, sep=' ')
            else:
                df.loc[subset].to_csv(path + name, header=None, mode='a', index=False, sep=' ')
                
    def tocsv(self, df, name, path=os.getcwd() + "/"):
        chunks = np.array_split(df.index, 100)
        for chunck, subset in enumerate(tqdm(chunks)):
            if chunck == 0: # first row
                df.loc[subset].to_csv(path + name, mode='w', index=True)
            else:
                df.loc[subset].to_csv(path + name, header=None, mode='a', index=True)

    # load DEG file from RNA-seq DEG analysis result
    def load_deg(self):
        try:
            print(">>>Loading DEG list file...")
            self.DEG = pd.read_csv("DEG_list.csv", sep=",")
            print(">>>DEG list file loaded...")
            return self.DEG
        except:
            print(">>>Fail to load DEG list file...")

    def load_gwas(self):  # load GWAS QC gwas file to acquire SNP position
        try:
            print(">>>Loading GWAS file...")
            self.gwas = pd.read_csv("BASEGWAS.assoc.logistic", header=0, sep="\t")
            self.gwas[self.gwas.columns[0]] = self.gwas[self.gwas.columns[0]].str.lstrip(" ")
            self.gwas = self.gwas[self.gwas.columns[0]].str.split(" +", expand=True)
            self.gwas.columns = ['CHR','SNP','BP','A1','TEST','NMISS','OR','SE','L95','U95','STAT','P']
            self.gwas = self.gwas[["CHR", "SNP", "BP"]]
            self.gwas.CHR = self.gwas.CHR.astype(str)
            print(">>>GWAS file loaded...")
            return self.gwas
        except (Exception,):
            print(">>>Fail to load GWAS file...")

    def load_refgene(self, url='https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz'):
        try:
            print(">>>Pulling refgene request...")
            request = requests.get(url, allow_redirects=True)
            file_name = url.split('/')[-1]
            print(">>>Writing refgene file...")
            open(file_name, "wb").write(request.content)
            self.refgene = pd.read_csv(file_name, compression='gzip', header=None, sep='\t', quotechar='"',
                                       on_bad_lines='skip')
            self.refgene[8] = self.refgene[8].apply(lambda x: x.split(";")[0].lstrip('gene_id ').strip('"'))
            print(">>>Keeping only transcript...")
            self.refgene = self.refgene[self.refgene[2]=="transcript"]
            self.refgene = self.refgene[[8,0,3,4]]
            self.refgene.drop_duplicates(inplace = True)
            self.refgene.columns = ["Gene","CHR","Start","End"]
            self.refgene.CHR = self.refgene.CHR.str.lstrip("chr")
            self.refgene.CHR = self.refgene.CHR.str.split('_').str[0]
            self.refgene.CHR = self.refgene.CHR.astype(str)
            os.remove(file_name)
            print(">>>Refgene file created...")
            return self.refgene
        except (Exception,):
            print(">>>Fail to create refgene file...")

    def load_tf(self):
        try:
            print(">>>Loading TF list...")
            self.tf_list = []
            with open("TF_list.txt") as f:
                for line in f.readlines():
                    self.tf_list.append(line.rsplit())
            self.tf_list = [x for item in self.tf_list for x in item]
            print(">>>TF list loaded...")
            return self.tf_list
        except (Exception,):
            print(">>>Fail to load TF list...")

    def get_matrix_id(self, tf):
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
        print(">>>Start JASPAR scraping...")
        i = len(self.tf_list)
        for tfs in self.tf_list:
            try:
                motif = self.get_matrix_id(tfs)
                if motif is None:
                    print(f'>>>Can not find JASPAR MOTIF ID for {tfs}!')
                else:
                    self.tfyesmotif_dict.setdefault(tfs, motif)
                    print(f'>>>The JASPAR MOTIF ID for {tfs} is {motif}!')
            except (Exception,):
                self.tfnomotif_list.append(tfs)
                print(f'>>>Can not find JASPAR MOTIF ID for {tfs}!')
            i -= 1
            if i != 0:
                print(f'>>>{i} TFs left!')
            else:
                print('>>>Scraping done!!')
        print(f'>>>We found {len(self.tfyesmotif_dict)} MOTIFs from our uploaded TFs by scraping JASPAR!')
        print(f'>>>{len(self.tfnomotif_list)} other TFs failed, listed below:')
        print(self.tfnomotif_list)
        return self.tfnomotif_list

    def paste_chip_url(self):
        print(">>>Pasting JASPAR url...")
        self.tfmotifURL_dict = {}
        for i in self.tfyesmotif_dict.keys():
            motif_id = self.tfyesmotif_dict[i]
            motif_url = 'https://jaspar.genereg.net/download/data/2022/bed/' + motif_id + '.bed'
            self.tfmotifURL_dict.setdefault(i, motif_url)
        print(">>>Done Pasting JASPAR url...")
        print(self.tfmotifURL_dict)
        return self.tfmotifURL_dict

    def load_chip_file(self, url):
        request = requests.get(url, allow_redirects=True)
        file_name = url.split('/')[-1]
        open(file_name, "wb").write(request.content)
        chip = pd.read_csv(file_name, header=None, sep='\t', quotechar='"', on_bad_lines='skip')
        return chip

    def liftover(self, chrom, pos):
        converter = ChainFile('hg19ToHg38.over.chain.gz', 'hg18', 'hg38')
        return converter[chrom][pos]

    def lift_beds(self):
        print(">>>Creating big chip df...")
        self.big_chip = pd.DataFrame()
        print(">>>Scraping chip files...")
        for i, j in self.tfmotifURL_dict.items():
            chip = self.load_chip_file(j)
            if len(chip.columns) < 6:
                os.remove(j.split('/')[-1])
                pass
            else:
                chip.columns = ["CHR", "Start", "End", "POS", "Peak", "Strand"]
                chip['Start'] = chip['Start'].astype('int')
                chip['End'] = chip['End'].astype('int')
                if chip['POS'].str.contains('hg19').any():
                    subprocess.run(["CrossMap.py", "bed", "hg19ToHg38.over.chain.gz", j.split('/')[-1], "new"])
                    chip2 = pd.read_csv("new", header=None, sep='\t', quotechar='"', on_bad_lines='skip')
                    chip2["ID"] = i
                    chip2.columns = ["CHR", "Start", "End", "POS", "Peak", "Strand", "ID"]
                    chip2 = chip2[["ID", "CHR", "Start", "End"]]
                    self.big_chip = self.big_chip.append(chip2)
                    os.remove(j.split('/')[-1])
                else:
                    chip["ID"] = i
                    chip = chip[["ID", "CHR", "Start", "End"]]
                    self.big_chip = self.big_chip.append(chip)
                    os.remove(j.split('/')[-1])
        os.remove("new")
        os.remove("new.unmap")
        print(">>>Done scraping chip files...")
        self.big_chip['CHR'] = self.big_chip['CHR'].str.strip('chr')
        self.big_chip['CHR'] = self.big_chip['CHR'].str.split('_').str[0]
        self.big_chip = self.big_chip.reset_index()
        print(">>>Writing big chip file...")
        self.totxt(self.big_chip, "BigChip.bed")
        print(self.big_chip.head(5))
        return self.big_chip

    def map_pos(self, chr, bp, n):
        mask = self.big_chip.loc[self.big_chip.CHR.values == chr].query("Start+1-{1} < {0} < End+1+{1}".format(bp, n))
        if len(mask) != 0:
          result = str(mask["ID"].tolist()).strip("[]")
          return result + ", " +  str(len(mask))
        else:
          return 0

    def parallelize_function_wide(self, df):
        df["ChIP_wide"] = df.apply(lambda r: self.map_pos(r["CHR"], r["BP"], 25), axis=1)
        return df

    def parallelize_function_narrow(self, df):
        df["ChIP_narrow"] = df.apply(lambda r: self.map_pos(r["CHR"], r["BP"], 10), axis=1)
        return df

    def parallelize_function_valid(self, df):
        df["ChIP_valid"] = df.apply(lambda r: self.map_pos(r["CHR"], r["BP"],0), axis=1)
        return df

    def parallelize_dataframe(self, df, func):
        num_processes = mp.cpu_count()
        print(f">>>Preparing parallel processing...\n>>>CPU count: {num_processes}")
        print(">>>Splitting gwas file...")
        num_processes = 20
        df_split = np.array_split(df, num_processes)
        print(">>>Parallel processing...")
        with mp.Pool(num_processes) as p:
          results = list(tqdm(p.imap(func, df_split), total=len(df_split)))
          new_df = pd.concat(results, axis=0, ignore_index=True)
        return new_df

    def mapping_tfbs(self):
        print(">>>Start mapping tfbs to gwas file...")
        # self.gwas['CHR'] = 'chr' + self.gwas['CHR'].astype(str)
        print(">>>Mapping with wide range...")
        self.gwas = self.parallelize_dataframe(self.gwas, self.parallelize_function_wide)
        self.totxt(self.gwas.SNP[self.gwas['ChIP_wide'] != 0], "SNP_wide_list.txt")
        gwas1 = self.gwas[self.gwas['ChIP_wide'] != 0]
        print(">>>Mapping with narrow range...")
        gwas1 = self.parallelize_dataframe(gwas1, self.parallelize_function_narrow)
        self.totxt(gwas1.SNP[gwas1['ChIP_narrow'] != 0], "SNP_narrow_list.txt")
        gwas2 = gwas1[gwas1['ChIP_narrow'] != 0]
        print(">>>Mapping with motif accurately...")
        gwas2 = self.parallelize_dataframe(gwas2, self.parallelize_function_valid)
        self.totxt(gwas2.SNP[gwas2['ChIP_valid'] != 0], "SNP_valid_list.txt")
        self.gwas_result = pd.merge(gwas1, gwas2, how="left", on=["CHR","SNP","BP","ChIP_wide","ChIP_narrow"])
        print(">>>Filling NA...")
        self.gwas_result = self.gwas_result.fillna(0)
        print(">>>Writing gwas_result.txt file...")
        self.totxt(self.gwas_result, "gwas_result.txt")
        print(">>>Mapping tfbs to gwas file done...")
        return self.gwas_result
    
    def map_ref(self, chr, bp):
        mask = self.refgene.loc[self.refgene.CHR.values == chr].query("Start < {0} < End".format(bp))
        try:
            if len(mask) != 0:
                a = mask["Gene"].iloc[0]
                b = mask["Start"].iloc[0]
                c = mask["End"].iloc[0]
                d = bp-mask["Start"].iloc[0]
                result = [a,b,c,d]
                return result
                # return str(mask["Gene"].tolist()).strip("([])")
            else:
                mask2 = self.refgene.loc[self.refgene.CHR.values == chr].query("{0} < Start".format(bp))
                mask2 = mask2.sort_values(by=['Start'])
                a2,b2,c2,d2 = mask2["Gene"].iloc[0], mask2["Start"].iloc[0], mask2["End"].iloc[0], bp - mask2["Start"].iloc[0]
                result2 = [a2,b2,c2,d2]
                return str(result2).strip("[]")
        except:
            return False
        
    def parallelize_function_mref(self, df):
        df["Gene"] = df.apply(lambda r: self.map_ref(r["CHR"], r["BP"]), axis=1)
        return df
    
    def gene_annotation(self):
        print(">>>Start gene annotation...")
        self.gwas_result.CHR = self.gwas_result.CHR.astype(str)
        self.gwas_result = self.parallelize_dataframe(self.gwas_result, self.parallelize_function_mref)
        self.gwas_result['Gene']  = self.gwas_result['Gene'].str.strip("[]")
        self.gwas_result[['Gene Name', 'Gene Start', 'Gene End', 'Distance']]  = self.gwas_result.Gene.str.split(",", expand = True)
        self.gwas_result = self.gwas_result.drop(columns=['Gene'])
        self.gwas_result['Gene Name']  = self.gwas_result["Gene Name"].str.strip("'")
        print(">>>Writing refgene annoted file...")
        self.totxt(self.gwas_result, "gwas_result_anotted.txt")
        print(">>>Done gene annotation...")
        return self.gwas_result
        
    def deg_annotation(self):
        print(">>>Start DEG annotation...")
        DEG_list = self.DEG["Gene symbol"].values.tolist()
        self.gwas_result["DEG"] = self.gwas_result["Gene Name"].isin(DEG_list)
        # pattern = '|'.join(self.DEG["Gene symbol"].astype(str).tolist())
        # self.gwas_result["DEG"] = self.gwas_result["Gene Name"].str.contains(pattern)
        print(">>>Writing deg annoted file...")
        self.totxt(self.gwas_result, "gwas_result_anotted_deg.txt")
        print(">>>Done DEG annotation...")
        return self.gwas_result
      
    def deg_upstream_filter(self):
        print(">>>Start filtering SNPs by position and DEG...")
        final_data = self.gwas_result.loc[self.gwas_result["ChIP_valid"] != 0]
        final_data = final_data.dropna(subset = ["Distance"])
        final_data["Distance"] = pd.to_numeric(final_data["Distance"], errors='coerce')
        final_data = final_data.loc[(final_data["Distance"] <= 2000) & (final_data["Distance"] >= -8000)]
        final_data = final_data.loc[final_data["DEG"] == True]
        final_data = final_data.iloc[:,[0,1,2,5,6,7,8,9,10]]
        print(">>>Done filtering, writing .csv file...")
        self.tocsv(final_data, "valid_deg_upstream_snp.csv")

if __name__ == '__main__':
    mapping_result = Rmcga()
    result_load = input(">>>YES to run mapping, NO to load result file! (Y/N): ")
    if result_load == "Y":
        mapping_result.load_gwas()
        chip_load = input(">>>YES to scrape JASPAR, NO to load big chip file! (Y/N): ")
        if chip_load == "Y":
            mapping_result.load_tf()
            mapping_result.scrape_jaspar()
            mapping_result.paste_chip_url()
            mapping_result.lift_beds()
        else:
            print(">>>Loading big chip file...")
            mapping_result.big_chip = pd.read_csv("BigChip.bed", header=0, sep=" ")
        mapping_result.mapping_tfbs()
    else:
        mapping_result.gwas_result = pd.read_csv("gwas_result.txt", header=0, sep=" ")
    mapping_result.load_refgene()
    mapping_result.load_deg()
    mapping_result.gene_annotation()
    mapping_result.deg_annotation()
    print(">>>Writing all mapping result csv...")
    mapping_result.tocsv(mapping_result.gwas_result, "gwas_result.csv")
    mapping_result.deg_upstream_filter()
    print(">>>All code done!")
    print("--- %s seconds ---" % (time.time() - start_time))
