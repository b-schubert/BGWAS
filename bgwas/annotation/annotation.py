"""
functions and objects to
annotate a genome sequence

:author: Benjamin Schubert
"""
import pandas as pd

from collections import namedtuple
from bisect import bisect_left
from intervaltree import IntervalTree
from Bio import SeqIO

# An element in of the genome annotation
GenomeEntry = namedtuple('GenomeEntry', ['type', 'name', 'locus', 'product', 'protein_id', 'strand', 'start', 'end'])


class GenomeAnnotation(object):
    """
    represents a genbank file
    and allows to efficiently annotate
    positions of interest
    """
    COLUMNS = ["type", "name", "locus", "product", 'protein_id', "strand", "start", "end"]

    def __init__(self, genbank_file):
        """
        initializes the GenomeAnnotation object

        :param genbank_file: a path to a genbank file
        """
        self.genome_tree = IntervalTree()
        self.gene_dic = {}
        self.locus_dic = {}
        self.type_dic = {}
        self.genome_id = None
        self.length = None
        self.__read_genbank(genbank_file)

        # internal data structure for quick internal nearest gene search if position is not annotated
        tmp = []
        for v in (self.type_dic["CDS"] + self.type_dic["gene"]):
            tmp.extend([(v.start, v), (v.end, v)])

        tmp.sort(key=lambda x: x[0])
        self.__index_list = []
        self.__cds_list = []
        for pos, cds in tmp:
            self.__index_list.append(pos)
            self.__cds_list.append(cds)

    def __read_genbank(self, genbank_file):
        """
        reads the genbank file and stores its content in a interval tree
        and other searchable containers for efficient querying

        :param genbank_file: a path to a genbank file
        """
        ##print("old implementation")
        pseudogenes = []
        with open(genbank_file, "r") as f:
            my_type, name, locus, product, product_id, strand, start, end = None, None, None, None, None, None, None, None

            annotated_features = set()

            # states
            gathering = False
            comment_block = False
            annotation_block = False
            c = 0
            for l in f:
                # skip empty lines
                if l.strip() == "":
                    continue

                splits = l.split()

                if splits[0].startswith("LOCUS"):
                    ##print(splits)
                    self.genome_id = splits[1].strip()
                    self.length = int(splits[2].strip())


                # are we at the end of the annotation block?
                if splits[0].startswith("ORIGIN"):
                    break

                # check for parsing stage
                if splits[0].startswith("COMMENT"):
                     comment_block = True

                if splits[0].startswith("FEATURES"):
                    ##print(annotated_features)
                    annotation_block = True
                    comment_block = False

                # COMMENT block feature annotation
                if comment_block and splits[0].startswith("Fe"):

                    gathering = True
                    for an in splits[3:]:
                        if not an.startswith("Gene"):
                            annotated_features.add(an.split(";")[0])
                        else:
                            annotated_features.add("gene")

                # FEATURES Block here we found an entry that we want to gather
                if annotation_block and splits[0] in annotated_features and ".." in splits[1]:

                    # first add already gathered entry into data structures
                    if locus is not None:
                        entry = GenomeEntry(my_type, name, locus, product, product_id, strand, start, end)

                        #if my_type == "PROMOTER":
                        #    print(entry)
                        # if its a gene annotation than first store it in temp for alter processing
                        if my_type == "gene":
                            pseudogenes.append(entry)
                        else:
                            if start > end:
                                ##print(entry)
                                c += 1
                                self.genome_tree.addi(start, self.length, entry)
                                self.genome_tree.addi(0, end, entry)
                            else:
                                self.genome_tree.addi(start, end, entry)

                            self.locus_dic[locus] = entry
                            self.type_dic.setdefault(my_type, []).append(entry)

                            if name is not None:
                                self.gene_dic[name] = entry

                        my_type, name, locus, product, product_id, strand, start, end = None, None, None, None, None, None, None, None

                    gathering = True
                    my_type = splits[0]
                    # determine strand, start and end

                    if splits[1].startswith('comp'):
                        interval = splits[1].strip('complement()')
                        strand = '-'
                    else:
                        interval = splits[1]
                        strand = '+'
                    start, end = map(lambda x: int(x) - 1, interval.split('..'))
                    # TODO: this has to be fixed in the genbank file
                    if start == end:
                        end += 1

                # gather annotated elements
                if gathering:

                    # if we are in the comment block than we are gathering annotated features
                    if comment_block:
                        if "::" in splits:
                            gathering = False
                        else:
                            for s in splits:
                                annotated_features.add(s.split(";")[0])

                    # if we are in the annotation block than we gather infos distributed over multiple lines
                    if annotation_block:
                        if splits[0].startswith("/locus"):
                            locus = l.split("=")[-1].replace('"', '').replace("_", "").strip()
                        elif splits[0].startswith("/product"):
                            product = l.split("=")[-1].replace('"', '').strip()
                        elif splits[0].startswith("/gene"):
                            name = l.split("=")[-1].replace('"', '').strip()
                        elif splits[0].startswith("/protein_id"):
                            product_id = l.split("=")[-1].replace('"', '').strip()
                        else:
                            continue

            # end of file
            if locus is not None:
                entry = GenomeEntry(my_type, name, locus, product, product_id, strand, start, end)
                # if its a gene annotation than first store it in temp for alter processing
                #if my_type == "PROMOTER":
                #    print(entry)
                if my_type == "gene":
                    pseudogenes.append(entry)
                else:
                    start = entry.start
                    end = entry.end
                    if start > end:
                        ##print(entry)
                        c +=1
                        self.genome_tree.addi(start, self.length, entry)
                        self.genome_tree.addi(0, end, entry)
                    else:
                        self.genome_tree.addi(entry.start, entry.end, entry)

                    self.locus_dic[locus] = entry
                    self.type_dic.setdefault(type, []).append(entry)

                    if name is not None:
                        self.gene_dic[name] = entry
            ##print("Wrongly start end", c)
            for p in pseudogenes:
                # if this is true gene did not have another entry
                if p.locus not in self.locus_dic:
                    self.locus_dic[p.locus] = p
                    self.type_dic.setdefault(p.type, []).append(p)
                    self.genome_tree.addi(p.start, p.end, p)
                    if p.name is not None:
                        self.gene_dic[p.name] = p

    def _read_genbank2(self, genbank_file):

        gene_tmp = []
        nop = [None]
        with open(genbank_file, "r") as gbk:
            anno = SeqIO.read(gbk, "genbank")
            self.genome_id = anno.id
            self.length = len(anno)

            for rec in anno.features:
                if rec.type == "source":
                    continue
                else:
                    entry = GenomeEntry(rec.type,
                                        rec.qualifiers.get("gene", nop)[0],
                                        rec.qualifiers.get("locus_tag", nop)[0],
                                        rec.qualifiers.get("product", nop)[0],
                                        rec.qualifiers.get("protein_id", nop)[0],
                                        "+" if rec.strand else "-",
                                        int(rec.location.start) - 1,
                                        int(rec.location.end) - 1)
                    if entry.type == "gene":
                        gene_tmp.append(entry)
                    else:
                        start = entry.start
                        end = entry.end
                        if start > end:
                            self.genome_tree.addi(start, self.length, entry)
                            self.genome_tree.addi(0, end, entry)
                        else:
                            self.genome_tree.addi(entry.start, entry.end, entry)

                        self.locus_dic[entry.locus] = entry
                        self.type_dic.setdefault(entry.type, []).append(entry)
                        if entry.name is not None:
                            self.gene_dic[entry.name] = entry

            for p in gene_tmp:
                # if this is true gene did not have another entry
                if p.locus not in self.locus_dic:
                    self.locus_dic[p.locus] = p
                    self.type_dic.setdefault(p.type, []).append(p)
                    self.genome_tree.addi(p.start, p.end, p)
                    if p.name is not None:
                        self.gene_dic[p.name] = p

    def __str__(self):
        return pd.DataFrame.from_records(list(self.locus_dic.values()), columns=self.COLUMNS).to_string()
    
    def annotate_positions(self, idx, aggregate=False):
        """
        annotates a list of positions with their associated genomic entries
        and returns a pandas dataframe with rows:

        pos, type, locus, name, product, strand, closest, distance, protein_pos, codon_pos

        :param idx: list of indices
        :return: pandas dataframe
        """

        # test if parameter is an iterable or int
        if isinstance(idx, int):
            idx = [idx]
        else:
            idx = list(set(idx))

        unknown = GenomeEntry("?", None, None, None, None, None, None, None)
        entries = []
        closest = []
        distance = []
        index = []
        
        protein_position = []
        codon_position = []
        
        for i in idx:
            data = self.genome_tree.search(i, strict=True)
            if data:
                # possible overlap of gene entries?
                for p in data:
                    #print(i, p.data)
                    index.append(i)
                    entries.append(p.data)
                    closest.append(None)
                    distance.append(None)
                    # calculate position within protein and codon position (1-indexed).
                    if p.data.strand == '+':
                        my_prot_pos = int((i - p.data.start) / 3) + 1 # int() rounds down.
                        my_codon_pos = ((i - p.data.start) % 3) + 1
                        ##print(my_prot_pos)
                    elif p.data.strand == '-':
                        my_prot_pos = int((p.data.end - i) / 3) + 1 # int() rounds down.
                        ##print(my_prot_pos)
                        my_codon_pos = ((p.data.end - i) % 3) + 1
                    else:
                        raise ValueError("strand annotation is invalid for gene, {}".format(p.data.locus))
                    protein_position.append(my_prot_pos)
                    codon_position.append(my_codon_pos)
            else:
                # position is not annotated in GenomeAnnotation
                # find closest annotated CDS
                index.append(i)
                entries.append(unknown)
                i_clos = self.find_closest_gene(i)
                closest.append(i_clos.locus)
                distance.append(min(abs(i - i_clos.start), abs(i - i_clos.end)))
                protein_position.append(None)
                codon_position.append(None)

        anno_df = pd.DataFrame.from_records(entries, columns=self.COLUMNS)

        anno_df["pos"] = index
        anno_df["closest"] = closest
        anno_df["distance"] = distance

        anno_df["protein_pos"] = protein_position
        anno_df["codon_pos"] = codon_position
        
        if aggregate:
            anno_df = anno_df.groupby("pos").agg(lambda col: ';'.join(map(str, col)))
            anno_df.reset_index(inplace=True)
            print(anno_df.head())
        return anno_df[["pos", "type", "locus", "name", "product",  "protein_id",
                        "strand", "closest", "distance", "start", "end", "protein_pos", "codon_pos"]]

    def find_closest_gene(self, pos):
        """
        Returns closest value to pos.
        If two numbers are equally close, return the smallest number.

        :param pos: the genome position
        :return: GenomeEntry
        """
        idx = bisect_left(self.__index_list, pos)
        if idx == 0:
            return self.__cds_list[0]
        if idx == len(self.__index_list):
            return self.__cds_list[-1]
        before = self.__index_list[idx - 1]
        after = self.__index_list[idx]
        if after - pos < pos - before:
            return self.__cds_list[idx]
        else:
            return self.__cds_list[idx - 1]

    def annotate_genes(self, genes):
        """
        annotates a list of gene and returns a pandas dataframe
        with the following columns:

        type name locus product strand start end

        :param genes: list of genes names
        :return: pandas dataframe
        """
        if isinstance(genes, str):
            genes = [genes]

        entries = [self.gene_dic[g] for g in genes if g in self.gene_dic]
        return pd.DataFrame.from_records(entries, columns=self.COLUMNS)

    def annotate_loci(self, loci):
        """
        annotates a list of loci tags and returns a pandas dataframe
        with the following columns:

        type name locus product strand start end

        :param loci: list of locus names
        :return: pandas dataframe
        """
        if isinstance(loci, str):
            loci = [loci]

        entries = [self.locus_dic[g] for g in loci if g in self.locus_dic]
        return pd.DataFrame.from_records(entries, columns=self.COLUMNS)

    def annotate_type(self, types):
        """
        annotates a list of types  and returns a pandas dataframe
        with the following columns:

        type name locus product strand start end

        :param types: list of types
        :return: pandas dataframe
        """
        if isinstance(types, str):
            types = [types]

        entries = []
        for g in types:
            if g in self.type_dic:
                for e in self.type_dic[g]:
                    entries.append(e)
        return pd.DataFrame.from_records(entries, columns=self.COLUMNS)

    def annotate_dataframe(self, df, column, suffix=("_x", "_y"), aggregate=False):
        """
        annotate an existing dataframe

        :param df: data frame to which annotation is added
        :param column: specifies the genome position column
        :param suffix: tuple of suffix that is added overlapping column names (default: (_x, _y))
        :param aggregate: determines whether duplicated entry are aggregated as a semicolon separated string
        :return: pandas dataframe
        """
        idx = set(df[column])
        pos_df = self.annotate_positions(idx, aggregate=aggregate)

        df = df.merge(pos_df, left_on=column, right_on="pos", how="inner", suffixes=suffix)
        df.drop("pos", axis=1, inplace=True)

        return df

if __name__ == "__main__":
    ## NOTE: I think Benni forgot to git add ../../test_data.gbk (or perhaps not).
    #g = GenomeAnnotation("../../test_data.gbk")
    g = GenomeAnnotation("/Users/Rohandinho/Dropbox (HMS)/Antibiotic Resistance/Analysis/N.G._full_b100k_k500/data/ng_FA1090_2016_with_promoters.gbk")
    #r = g.annotate_positions(160)
    r = g.annotate_positions(1668203)
    r.groupby("pos")
    print(r.groupby("pos").agg(lambda col: ';'.join(map(str, col))))

