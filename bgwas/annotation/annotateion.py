"""
functions and objects to
annotate a genome sequence

:author: Benjamin Schubert
"""
import pandas as pd

from collections import namedtuple
from bisect import bisect_left
from intervaltree import IntervalTree


# An element in of the genome annotation
GenomeEntry = namedtuple('GenomeEntry', ['type', 'name', 'locus', 'product', 'strand', 'start', 'end'])


class GenomeAnnotation(object):
    """
    represents a genbank file
    and allows to efficiently annotate
    positions of interest
    """
    COLUMNS = ["type", "name", "locus", "product", "strand", "start", "end"]

    def __init__(self, genbank_file):
        """
        initializes the GenomeAnnotation object

        :param genbank_file: a path to a genbank file
        """
        self.genome_tree = IntervalTree()
        self.gene_dic = {}
        self.locus_dic = {}
        self.type_dic = {}

        self.__read_genbank(genbank_file)

        # internal data structure for quick internal nearest gene search if position is not annotated
        tmp = []

        for v in self.type_dic["CDS"]:
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

        with open(genbank_file, "r") as f:
            type, name, locus, product, strand, start, end = None, None, None, None, None, None, None

            annotated_features = set()

            # states
            gathering = False
            comment_block = False
            annotation_block = False

            for l in f:
                # skip empty lines
                if l.strip() == "":
                    continue

                splits = l.split()

                # are we at the end of the annotation block?
                if splits[0].startswith("ORIGIN"):
                    break

                # check for parsing stage
                if splits[0].startswith("COMMENT"):
                    comment_block = True

                if splits[0].startswith("FEATURES"):
                    annotation_block = True
                    comment_block = False

                # COMMENT block feature annotation
                if comment_block and splits[0].startswith("Fe"):
                    gathering = True
                    for an in splits[3:]:
                        # we don't want to gather the gene annotation - redundant to CDS and some times contains
                        # pseudo genes
                        if not an.startswith("Gene"):
                            annotated_features.add(an.split(";")[0])

                # FEATURES Block here we found an entry that we want to gather
                if annotation_block and splits[0] in annotated_features and ".." in splits[1]:

                    # first add already annotated entry into data structures
                    if locus is not None:
                        entry = GenomeEntry(type, name, locus, product, strand, start, end)
                        self.locus_dic[locus] = entry
                        self.type_dic.setdefault(type, []).append(entry)
                        self.genome_tree.addi(start, end, entry)
                        if name is not None:
                            self.gene_dic[name] = entry

                        type, name, locus, product, strand, start, end = None, None, None, None, None, None, None

                    gathering = True
                    type = splits[0]
                    # determine strand, start and end
                    if splits[1].startswith('comp'):
                        interval = splits[1].strip('complement()')
                        strand = '-'
                    else:
                        interval = splits[1]
                        strand = '+'
                    try:
                        start, end = map(lambda x: int(x) - 1, interval.split('..'))
                    except:
                        print(splits)
                        print(splits[1], interval)
                        break

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
                        else:
                            continue

            # end of file
            if locus is not None:
                entry = GenomeEntry(type, name, locus, product, strand, start, end)
                self.locus_dic[locus] = entry
                self.type_dic.setdefault(type, []).append(entry)
                self.genome_tree.addi(start, end, entry)
                if name is not None:
                    self.gene_dic[name] = entry

    def __str__(self):
        return pd.DataFrame.from_records(list(self.locus_dic.values()), columns=self.COLUMNS).to_string()

    def annotate_positions(self, idx):
        """
        annotates a list of positions with their associated genomic entries
        and returns a pandas dataframe with rows:

        pos, type, locus, name, product, strand, closest, distance

        :param idx: list of indices
        :return: pandas dataframe
        """

        # test if parameter is an iterable or int
        if isinstance(idx, int):
            idx = [idx]

        unknown = GenomeEntry("?", None, None, None, None, None, None)
        entries = []
        closest = []
        distance = []
        index = []
        for i in idx:
            data = self.genome_tree.search(i, strict=True)
            if data:
                # possible overlap of gene entries?
                p = data.pop()
                index.append(i)
                entries.append(p.data)
                closest.append(None)
                distance.append(None)
            else:
                # position is not annotated in GenomeAnnotation
                # find closest annotated CDS
                index.append(i)
                entries.append(unknown)
                i_clos = self.find_closest_gene(i)
                closest.append(i_clos.locus)
                distance.append(min(abs(i - i_clos.start), abs(i - i_clos.end)))
        print("entry", len(entries), "idx", len(index), "closest", len(closest))
        anno_df = pd.DataFrame.from_records(entries, columns=self.COLUMNS)

        anno_df["pos"] = index
        anno_df["closest"] = closest
        anno_df["distance"] = distance

        return anno_df[["pos", "type", "locus", "name", "product", "strand", "closest", "distance"]]

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

    def annotate_dataframe(self, df, column, suffix=("_x", "_y")):
        """
        annotate an existing dataframe

        :param df: data frame to which annotation is added
        :param column: specifies the genome position column
        :param suffix: tuple of suffix that is added overlapping column names (default: (_x, _y))
        :return: pandas dataframe
        """
        idx = set(df[column])
        pos_df = self.annotate_positions(idx)

        df = df.merge(pos_df, left_on=column, right_on="pos", how="inner", suffixes=suffix)
        df.drop("pos", axis=1, inplace=True)

        return df


