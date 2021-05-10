# Loading necessary modules
# module for retrieving directory name
import os

# modules for data manipulation
import pandas as pd
import numpy as np
import vcf
import json
import logging

# module for command line prompt
import argparse

# modules from the plotting library bokeh
from bokeh.layouts import row, column, layout
from bokeh.models import Span, Label, OpenURL, TapTool, NumberFormatter, CustomJS, MultiSelect, DataRange1d
from bokeh.models.widgets import DataTable, TableColumn, Div
from bokeh.plotting import *
from bokeh.transform import factor_cmap
from bokeh.resources import INLINE

parser = argparse.ArgumentParser(description="Visualize CNV data from short read sequencing data.")

parser.add_argument("--ratio-file", "-r", required=True,
                    dest="ratio_file", default=None,
                    help="File which contains the log2(FC) of bins between the tumor sample and another normal sample."
                         " [Required]")

parser.add_argument("--genome-file", "-x", required=True,
                    dest="genome_file", default=None,
                    help="File which contains chromosome length and cumulative genomic length."
                         " [Required]")

parser.add_argument("--config-file", "-c", required=True,
                    dest="config_file", default=None,
                    help="File which contains plot options and column name customizations."
                         " [Required]")

parser.add_argument("--out-dir", "-d", required=True,
                    dest="out_dir", default=None,
                    help="Directory to place output files."
                         " [Required]")

parser.add_argument("--out-file", "-o", required=True,
                    dest="out_file", default=None,
                    help="Output file name (file will be placed in output dir - enter only filename)."
                         " [Required]")

parser.add_argument("--seg-file", "-s", required=False,
                    dest="seg_file", default=None,
                    help="File which contains the segmentation of log2(FC) bin values between the tumor sample "
                         "and normal sample.")

parser.add_argument("--gene-file", "-g", required=False,
                    dest="gene_file", default=None,
                    help="File which contains gene calling information.")

parser.add_argument("--seg-blacklist", "-t", required=False,
                    dest="seg_blacklist", default=None,
                    help="BED file of problematic copy number regions to highlight.")

parser.add_argument("--annotation-file", "-a", required=False,
                    dest="annot_file", default=None,
                    help="File which contains gene/exon information.")

parser.add_argument("--vcf-file", "-v", required=False,
                    dest="vcf_file", default=None,
                    help="VCF containing variants to plot VAF.")

parser.add_argument("--recenter", "-y", required=False,
                    dest="recenter", default=None,
                    help="Recenter to provided log2(FC).")

parser.add_argument("--vcf-filt-file", "-f", required=False,
                    dest="vcf_filt_file", default=None, action="store_true",
                    help="Flag to output filtered variants used for plotting VAFs. (applicable only if providing VCF)")

parser.add_argument("--vcf-blacklist", "-b", required=False,
                    dest="bed_blacklist", default=None,
                    help="File containing variants to NOT plot VAF. (applicable only if providing VCF)")

parser.add_argument("--purity", "-p", required=False,
                    dest="purity", default=None,
                    help="Purity of the sample.")

parser.add_argument("--ploidy", "-l", required=False,
                    dest="ploidy", default=None,
                    help="Ploidy of the sample.")

parser.add_argument("--gender", "-z", required=False,
                    dest="gender", default=None,
                    help="Ploidy of the sample.")

parser.add_argument("--verbose", "-j", required=False,
                    dest="verb_log", default=None, action="store_true",
                    help="Verbose logging output")

parser.add_argument('--version', action='version', version='%(prog)s v1.0.0')

options = parser.parse_args()

if options.verb_log:
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.DEBUG)
else:
    logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

outdir = options.out_dir

with open(options.config_file, "r") as json_file:
    config = json.load(json_file)

logging.info("Successfully read the configuration file.")


def draw_chr_boundary(figure, chr_boundary, genome, vaf):
    # first line drawn at the beginning of the plotting panel
    chr_temp_line = Span(location=0,
                         dimension='height',
                         line_color=config['plots']['chromosome_boundaries']['line_color'],
                         line_width=config['plots']['chromosome_boundaries']['line_width'],
                         line_dash=config['plots']['chromosome_boundaries']['line_dash'],
                         line_alpha=config['plots']['chromosome_boundaries']['line_alpha'])
    figure.add_layout(chr_temp_line)
    temp_loc = 0

    # check if to be drawn by the genome coordinates or by data points
    if (genome == False):
        chr_boundary = chr_boundary.sort_values(by="ind")
        # for loop to draw subsequent chromosome boundary lines
        for index, row in chr_boundary.iterrows():
            temp_text_loc = round(temp_loc + (row['ind'] - temp_loc) / 2)
            if (vaf == True):
                text = Label(x=temp_text_loc, y=0.1,
                             text=str(row[config['files']['ratio_file']['column_names']['chromosome']]),
                             text_color=config['plots']['chromosome_boundaries']['text_color'],
                             text_font_size=config['plots']['chromosome_boundaries']['text_font_size'],
                             background_fill_alpha=config['plots']['chromosome_boundaries']['text_background_alpha'],
                             background_fill_color=config['plots']['chromosome_boundaries']['text_background_color'],
                             text_font_style=config['plots']['chromosome_boundaries']['text_font_style'])
            else:
                text = Label(x=temp_text_loc, y=-2.5,
                             text=str(row[config['files']['ratio_file']['column_names']['chromosome']]),
                             text_color=config['plots']['chromosome_boundaries']['text_color'],
                             text_font_size=config['plots']['chromosome_boundaries']['text_font_size'],
                             background_fill_alpha=config['plots']['chromosome_boundaries']['text_background_alpha'],
                             background_fill_color=config['plots']['chromosome_boundaries']['text_background_color'],
                             text_font_style=config['plots']['chromosome_boundaries']['text_font_style'])
            figure.add_layout(text)

            chr_temp_line = Span(location=row['ind'],
                                 dimension='height',
                                 line_color=config['plots']['chromosome_boundaries']['line_color'],
                                 line_width=config['plots']['chromosome_boundaries']['line_width'],
                                 line_dash=config['plots']['chromosome_boundaries']['line_dash'],
                                 line_alpha=config['plots']['chromosome_boundaries']['line_alpha'])
            figure.add_layout(chr_temp_line)
            temp_loc = row['ind']
    else:
        chr_boundary = chr_boundary.sort_values(by="genome_cumsum")
        # for loop to draw subsequent chromosome boundary lines
        for index, row in chr_boundary.iterrows():
            temp_text_loc = round(temp_loc + (row['genome_cumsum'] - temp_loc) / 2)
            if (vaf == True):
                text = Label(x=temp_text_loc, y=0.1,
                             text=str(row[config['files']['ratio_file']['column_names']['chromosome']]),
                             text_color=config['plots']['chromosome_boundaries']['text_color'],
                             text_font_size=config['plots']['chromosome_boundaries']['text_font_size'],
                             background_fill_alpha=config['plots']['chromosome_boundaries']['text_background_alpha'],
                             background_fill_color=config['plots']['chromosome_boundaries']['text_background_color'],
                             text_font_style=config['plots']['chromosome_boundaries']['text_font_style'])
            else:
                text = Label(x=temp_text_loc, y=-2.5,
                             text=str(row[config['files']['ratio_file']['column_names']['chromosome']]),
                             text_color=config['plots']['chromosome_boundaries']['text_color'],
                             text_font_size=config['plots']['chromosome_boundaries']['text_font_size'],
                             background_fill_alpha=config['plots']['chromosome_boundaries']['text_background_alpha'],
                             background_fill_color=config['plots']['chromosome_boundaries']['text_background_color'],
                             text_font_style=config['plots']['chromosome_boundaries']['text_font_style'])
            figure.add_layout(text)

            chr_temp_line = Span(location=row['genome_cumsum'],
                                 dimension='height',
                                 line_color=config['plots']['chromosome_boundaries']['line_color'],
                                 line_width=config['plots']['chromosome_boundaries']['line_width'],
                                 line_dash=config['plots']['chromosome_boundaries']['line_dash'],
                                 line_alpha=config['plots']['chromosome_boundaries']['line_alpha'])
            figure.add_layout(chr_temp_line)
            temp_loc = row['genome_cumsum']


# read ratio file for plotting log2(FC) points
data = pd.read_csv(options.ratio_file, sep=config['files']['ratio_file']['field_separator'])
logging.info("Successfully read the ratio file.")

# make sure chromosome names are string data type
data[config['files']['ratio_file']['column_names']['chromosome']] = \
    data[config['files']['ratio_file']['column_names']['chromosome']].astype(str)

# index generation for marking chromosomes and drawing chromosome lines
data["ind"] = range(len(data))

if (options.recenter):
    data[config['files']['ratio_file']['column_names']['log2FC']] = \
        data[config['files']['ratio_file']['column_names']['log2FC']] - float(options.recenter)

# check if gene annotation column is provided accurately
if config['files']['ratio_file']['column_names']['gene'] not in data.columns:
    data[config['files']['ratio_file']['column_names']['gene']] = "-"
    logging.warning("Column titled \"" + config['files']['ratio_file']['column_names'][
        'gene'] + "\" not found in the ratio file. Using \"-\" in place of the annotation. "
                  "Even off-target bins might be marked as on-target!")

# check if weight column is provided accurately
if config['files']['ratio_file']['column_names']['weight'] not in data.columns:
    data[config['files']['ratio_file']['column_names']['weight']] = 0.1
    logging.warning("Column titled \"" + config['files']['ratio_file']['column_names'][
        'weight'] + "\" not found in the ratio file. Using \"0.1\" in place of the weight.")

# color background off-target bins (aka Antitarget) differently from on-target bins (aka Target)
data['label'] = np.where(
    (data[config['files']['ratio_file']['column_names']['gene']] == config['files']['ratio_file']['off_target_label']),
    "Antitarget", "Target")

# do not plot antitarget bins that are below log2 -10 (low confidence points)
data[config['files']['ratio_file']['column_names']['log2FC']] = np.where(np.logical_and(
    data[config['files']['ratio_file']['column_names']['log2FC']] < config['files']['ratio_file'][
        'off_target_low_conf_log2'], data.label == "Antitarget"), np.nan, data[config['files']['ratio_file'][
    'column_names']['log2FC']])

if data.label[0] == "Antitarget" and np.isnan(data[config['files']['ratio_file']['column_names']['log2FC']][0]):
    data.loc[0, config['files']['ratio_file']['column_names']['log2FC']] = 0
    #print("Changing first log2FC value manually to 0! -> related to Bokeh not able to serialize beginning with NaN")

chr_cumsum = pd.read_csv(options.genome_file, sep=config['files']['genome_file']['field_separator'])
logging.info("Successfully read the genome file.")

# make sure chromosome names are string data type
chr_cumsum[config['files']['genome_file']['column_names']['chromosome']] = \
    chr_cumsum[config['files']['genome_file']['column_names']['chromosome']].astype(str)

# retrieve genome level coordinates from chromosome level coordinates
data = pd.merge(data, chr_cumsum, left_on=config['files']['ratio_file']['column_names']['chromosome'],
                right_on=config['files']['genome_file']['column_names']['chromosome'], how='left')
data['genome_cumsum'] = data[config['files']['ratio_file']['column_names']['start']] + data[
    config['files']['genome_file']['column_names']['chr_cumulative_length']]

if (options.seg_blacklist):
    seg_blacklist = pd.read_csv(options.seg_blacklist, sep="\t", header=None)
    logging.info("Successfully read the BED segment artifact file.")

    seg_blacklist = seg_blacklist[seg_blacklist.columns[0:5]]
    seg_blacklist.columns = ['chromosome', 'start', 'end']
    seg_blacklist['chromosome'] = seg_blacklist['chromosome'].astype(str)

    seg_blacklist = pd.merge(seg_blacklist, chr_cumsum, left_on="chromosome",
                   right_on=config['files']['genome_file']['column_names']['chromosome'], how='left')
    seg_blacklist['genome_cumsum_start'] = seg_blacklist["start"] + seg_blacklist[
        config['files']['genome_file']['column_names']['chr_cumulative_length']]
    seg_blacklist['genome_cumsum_end'] = seg_blacklist["end"] + seg_blacklist[
        config['files']['genome_file']['column_names']['chr_cumulative_length']]
    seg_blacklist['index'] = range(len(seg_blacklist))

if (options.seg_file):
    seg = pd.read_csv(options.seg_file, sep=config['files']['segmentation_file']['field_separator'])
    logging.info("Successfully read the genome file.")

    if (options.recenter):
        seg[config['files']['segmentation_file']['column_names']['log2FC']] = \
            seg[config['files']['segmentation_file']['column_names']['log2FC']] - float(options.recenter)

    # check if gene annotation column is provided accurately
    if config['files']['segmentation_file']['column_names']['gene'] not in seg.columns:
        seg[config['files']['segmentation_file']['column_names']['gene']] = "-"
        logging.warning("Column titled \"" + config['files']['segmentation_file']['column_names'][
            'gene'] + "\" not found in the segmentation file. Using \"-\" in place of the annotation.")

    # make sure chromosome names are string data type
    seg[config['files']['segmentation_file']['column_names']['chromosome']] = \
        seg[config['files']['segmentation_file']['column_names']['chromosome']].astype(str)
    seg["prob_start"] = data.merge(seg, how="right",
                                   left_on=[config['files']['ratio_file']['column_names']['chromosome'],
                                            config['files']['ratio_file']['column_names']['start']],
                                   right_on=[config['files']['segmentation_file']['column_names']['chromosome'],
                                             config['files']['segmentation_file']['column_names']['start']])[
        "ind"]
    seg["prob_end"] = data.merge(seg, how="right", left_on=[config['files']['ratio_file']['column_names']['chromosome'],
                                                            config['files']['ratio_file']['column_names']['end']],
                                 right_on=[config['files']['segmentation_file']['column_names']['chromosome'],
                                           config['files']['segmentation_file']['column_names']['end']])[
        "ind"]
    seg = pd.merge(seg, chr_cumsum, left_on=config['files']['segmentation_file']['column_names']['chromosome'],
                   right_on=config['files']['genome_file']['column_names']['chromosome'], how='left')
    seg['genome_cumsum_start'] = seg[config['files']['segmentation_file']['column_names']['start']] + seg[
        config['files']['genome_file']['column_names']['chr_cumulative_length']]
    seg['genome_cumsum_end'] = seg[config['files']['segmentation_file']['column_names']['end']] + seg[
        config['files']['genome_file']['column_names']['chr_cumulative_length']]
    seg['index'] = range(len(seg))

if (options.gene_file):
    gen = pd.read_csv(options.gene_file, sep=config['files']['gene_file']['field_separator'])
    logging.info("Successfully read the gene CNV file.")

    if (options.recenter):
        gen[config['files']['gene_file']['column_names']['log2FC']] = \
            gen[config['files']['gene_file']['column_names']['log2FC']] - float(options.recenter)

    gen[config['files']['gene_file']['column_names']['chromosome']] = \
        gen[config['files']['gene_file']['column_names']['chromosome']].astype(str)
    gen["prob_start"] = data.merge(gen, how="right",
                                   left_on=[config['files']['ratio_file']['column_names']['chromosome'],
                                            config['files']['ratio_file']['column_names']['start']],
                                   right_on=[config['files']['gene_file']['column_names']['chromosome'],
                                             config['files']['gene_file']['column_names']['start']])[
        "ind"]
    gen["prob_end"] = data.merge(gen, how="right", left_on=[config['files']['ratio_file']['column_names']['chromosome'],
                                                            config['files']['ratio_file']['column_names']['end']],
                                 right_on=[config['files']['gene_file']['column_names']['chromosome'],
                                           config['files']['gene_file']['column_names']['end']])["ind"]
    gen = pd.merge(gen, chr_cumsum, left_on=config['files']['gene_file']['column_names']['chromosome'],
                   right_on=config['files']['genome_file']['column_names']['chromosome'], how='left')
    gen['genome_cumsum_start'] = gen[config['files']['gene_file']['column_names']['start']] + gen[
        config['files']['genome_file']['column_names']['chr_cumulative_length']]
    gen['genome_cumsum_end'] = gen[config['files']['gene_file']['column_names']['end']] + gen[
        config['files']['genome_file']['column_names']['chr_cumulative_length']]

    #subsetting columns to define gene boundaries
    gene_boundaries = gen[[config['files']['gene_file']['column_names']['chromosome'],
                           config['files']['gene_file']['column_names']['start'],
                           config['files']['gene_file']['column_names']['end'],
                           config['files']['gene_file']['column_names']['gene']]]
    gene_boundaries.columns = ['chromosome', 'gene_start', 'gene_end', 'gene']
    gene_boundaries = pd.merge(gene_boundaries, chr_cumsum, left_on='chromosome',
                               right_on=config['files']['genome_file']['column_names']['chromosome'], how='left')
    gene_boundaries['genome_cumsum_start'] = gene_boundaries.gene_start + gene_boundaries[
        config['files']['genome_file']['column_names']['chr_cumulative_length']]
    gene_boundaries['genome_cumsum_end'] = gene_boundaries.gene_end + gene_boundaries[
        config['files']['genome_file']['column_names']['chr_cumulative_length']]

if (options.annot_file):
    gene_track = pd.read_csv(options.annot_file, sep=config['files']['annotation_file']['field_separator'])
    logging.info("Successfully read the annotation file.")

    gene_track[config['files']['annotation_file']['column_names']['chromosome']] = \
        gene_track[config['files']['annotation_file']['column_names']['chromosome']].astype(str)
    gene_track = pd.merge(gene_track, chr_cumsum,
                          left_on=config['files']['annotation_file']['column_names']['chromosome'],
                          right_on=config['files']['genome_file']['column_names']['chromosome'], how='left')
    gene_track['genome_cumsum_start'] = gene_track[config['files']['annotation_file']['column_names']['exon_start']] + \
                                        gene_track[
                                            config['files']['genome_file']['column_names']['chr_cumulative_length']]
    gene_track['genome_cumsum_end'] = gene_track[config['files']['annotation_file']['column_names']['exon_end']] + \
                                      gene_track[
                                          config['files']['genome_file']['column_names']['chr_cumulative_length']]

if (options.bed_blacklist and options.vcf_file):
    bed_blacklist = pd.read_csv(options.bed_blacklist, sep="\t", header=None)
    logging.info("Successfully read the BED SNP blacklist file.")

    bed_blacklist = bed_blacklist[bed_blacklist.columns[0:5]]
    bed_blacklist.columns = ['chromosome', 'start', 'ref', 'alt', 'count']
    bed_blacklist['chromosome'] = bed_blacklist['chromosome'].astype(str)

if (options.vcf_file):
    df_vaf = []
    logging.info("Reading VCF file ...")
    vcf_reader = vcf.Reader(open(options.vcf_file, 'r'))
    logging.info("Processing VCF file ...")
    for record in vcf_reader:
        if record.INFO[config['files']['vcf_file']['info_fields']['depth']] >= \
                config['files']['vcf_file']['thresholds']['depth'] and \
                record.CHROM != 'X' and \
                len(record.REF) == 1 and \
                len(record.ALT) == 1 and \
                len(record.ALT[0]) == 1 and \
                ((record.INFO[config['files']['vcf_file']['info_fields']['forward_alt_reads']][0] +
                  record.INFO[config['files']['vcf_file']['info_fields']['reverse_alt_reads']][0]) / record.INFO[
                     config['files']['vcf_file']['info_fields']['depth']] >= config['files']['vcf_file']['thresholds'][
                     'low_vaf_filter'] and
                 (record.INFO[config['files']['vcf_file']['info_fields']['forward_alt_reads']][0] +
                  record.INFO[config['files']['vcf_file']['info_fields']['reverse_alt_reads']][0]) / record.INFO[
                     config['files']['vcf_file']['info_fields']['depth']] <= config['files']['vcf_file']['thresholds'][
                     'high_vaf_filter']) and \
                record.INFO[config['files']['vcf_file']['info_fields']['forward_alt_reads']][0] >= \
                config['files']['vcf_file']['thresholds']['forward_alt_reads'] and \
                record.INFO[config['files']['vcf_file']['info_fields']['reverse_alt_reads']][0] >= \
                config['files']['vcf_file']['thresholds']['reverse_alt_reads']:
            df_vaf.append({'chromosome': record.CHROM,
                           'start': record.POS,
                           'ref': record.REF,
                           'alt': record.ALT[0],
                           'DP': record.INFO[config['files']['vcf_file']['info_fields']['depth']],
                           'AF': (record.INFO[config['files']['vcf_file']['info_fields']['forward_alt_reads']][0] +
                                  record.INFO[config['files']['vcf_file']['info_fields']['reverse_alt_reads']][0]) /
                                 record.INFO[config['files']['vcf_file']['info_fields']['depth']]})
    df_vaf = pd.DataFrame(df_vaf)
    logging.info("Successfully filtered VCF file.")
    if (not df_vaf.empty):
        df_vaf['chromosome'] = df_vaf['chromosome'].astype(str)
        df_vaf['ref'] = df_vaf['ref'].astype('str')
        df_vaf['alt'] = df_vaf['alt'].astype('str')

        if (options.vcf_filt_file):
            output_VAF_filename = os.path.basename(os.path.splitext(options.vcf_file)[0]) + "_filt_SNPs.txt"
            print("Writing plotted variants to ", outdir, "/", output_VAF_filename, "...", sep=" ")
            df_vaf.to_csv(outdir + "/" + output_VAF_filename, sep="\t", index=False)

        # print(len(df_vaf))
        df_vaf = pd.merge(df_vaf, chr_cumsum, left_on='chromosome',
                          right_on=config['files']['genome_file']['column_names']['chromosome'], how='left')

        if (options.bed_blacklist):
            df_vaf = pd.merge(df_vaf, bed_blacklist, left_on=["chromosome", "start", "ref", "alt"],
                              right_on=["chromosome", "start", "ref", "alt"], how="left")
            df_vaf['count'].replace('', np.nan, inplace=True)
            df_vaf = df_vaf[pd.isnull(df_vaf['count'])]
        df_vaf['genome_cumsum'] = df_vaf.start + df_vaf[
            config['files']['genome_file']['column_names']['chr_cumulative_length']]

        if (options.seg_file):
            df_vaf_seg = (pd.merge(df_vaf, seg, left_on='chromosome',
                                   right_on=config['files']['segmentation_file']['column_names']['chromosome'],
                                   how='left'))
            if (config['files']['segmentation_file']['column_names']['start'] == "start"):
                df_vaf_seg['intersect'] = np.where(
                    ((df_vaf_seg['start_x'] >= df_vaf_seg['start_y']) & (df_vaf_seg['start_x'] <= df_vaf_seg[
                        config['files']['segmentation_file']['column_names']['end']])), 1, 0)
                df_vaf_seg = df_vaf_seg.query('intersect == 1')
                df_vaf_seg['cluster'] = np.where(df_vaf_seg['AF'] > 0.5, 1, 0)
                df_vaf_count = df_vaf_seg.groupby(
                    ['chromosome', 'start_y', 'end', 'cluster']).count().reset_index().intersect
                df_vaf_seg = df_vaf_seg.groupby(['chromosome', 'start_y', 'end', 'cluster']).mean().reset_index()
                df_vaf_seg['count'] = df_vaf_count
            else:
                df_vaf_seg['intersect'] = np.where(
                    ((df_vaf_seg['start'] >= df_vaf_seg[
                        config['files']['segmentation_file']['column_names']['start']]) & (
                                 df_vaf_seg['start'] <= df_vaf_seg[
                             config['files']['segmentation_file']['column_names']['end']])), 1, 0)
                df_vaf_seg = df_vaf_seg.query('intersect == 1')
                df_vaf_seg['cluster'] = np.where(df_vaf_seg['AF'] > 0.5, 1, 0)
                df_vaf_count = df_vaf_seg.groupby(
                    ['chromosome', config['files']['segmentation_file']['column_names']['start'],
                     config['files']['segmentation_file']['column_names']['end'],
                     'cluster']).count().reset_index().intersect
                df_vaf_seg = df_vaf_seg.groupby(
                    ['chromosome', config['files']['segmentation_file']['column_names']['start'],
                     config['files']['segmentation_file']['column_names']['end'], 'cluster']).mean().reset_index()
                df_vaf_seg['count'] = df_vaf_count

        if (options.gene_file):
            df_vaf_gene = (pd.merge(df_vaf, gene_boundaries, on='chromosome', how='left'))
            df_vaf_gene['intersect'] = np.where(
                ((df_vaf_gene['start'] >= df_vaf_gene['gene_start']) & (
                            df_vaf_gene['start'] <= df_vaf_gene['gene_end'])), 1, 0)
            df_vaf_gene = df_vaf_gene.query('intersect == 1')
            df_vaf_gene['cluster'] = np.where(df_vaf_gene['AF'] > 0.5, 1, 0)
            df_vaf_gene = df_vaf_gene.groupby(['chromosome', 'gene_start', 'gene_end', 'cluster', 'genome_cumsum_start',
                                               'genome_cumsum_end']).mean().reset_index()
    else:
        logging.warning("No variants to plot! (VCF did not have any variants passing the set thresholds)")

# setting up output HTML file
logging.info("Setting up output file.")
output_filename = os.path.basename(options.out_file)
title_name = os.path.basename(os.path.splitext(options.ratio_file)[0])
logging.info("Writing visualization to " + outdir + "/" + output_filename + "...")
output_file(outdir + "/" + output_filename, title=title_name + ' reconCNV Plot')

source = ColumnDataSource(data=dict(
    chrom=data[config['files']['ratio_file']['column_names']['chromosome']],
    start=data[config['files']['ratio_file']['column_names']['start']],
    genome_cumsum=data.genome_cumsum,
    end=data[config['files']['ratio_file']['column_names']['end']],
    logFC=data[config['files']['ratio_file']['column_names']['log2FC']],
    gene=data[config['files']['ratio_file']['column_names']['gene']],
    ind=data.ind,
    label=data.label,
    weight=data[config['files']['ratio_file']['column_names']['weight']] * config['files']['ratio_file'][
        'weight_scaling_factor']))

source_copy = ColumnDataSource(data=dict(
    chrom=data[config['files']['ratio_file']['column_names']['chromosome']],
    start=data[config['files']['ratio_file']['column_names']['start']],
    genome_cumsum=data.genome_cumsum,
    end=data[config['files']['ratio_file']['column_names']['end']],
    logFC=data[config['files']['ratio_file']['column_names']['log2FC']],
    gene=data[config['files']['ratio_file']['column_names']['gene']],
    ind=data.ind,
    label=data.label,
    weight=data[config['files']['ratio_file']['column_names']['weight']] * config['files']['ratio_file'][
        'weight_scaling_factor']))

# if both integer copy number and clonality information is present
if (options.seg_file and
        ((config['files']['segmentation_file']['column_names']['major_cn'] in seg.columns and
        config['files']['segmentation_file']['column_names']['minor_cn'] in seg.columns) and
        config['files']['segmentation_file']['column_names']['cell_frac'] in seg.columns)):
    source_seg = ColumnDataSource(data=dict(
        chrom=seg[config['files']['segmentation_file']['column_names']['chromosome']],
        start=seg[config['files']['segmentation_file']['column_names']['start']],
        end=seg[config['files']['segmentation_file']['column_names']['end']],
        logFC=seg[config['files']['segmentation_file']['column_names']['log2FC']],
        prob_start=seg.prob_start,
        prob_end=seg.prob_end,
        genome_cumsum_start=seg.genome_cumsum_start,
        genome_cumsum_end=seg.genome_cumsum_end,
        major_cn=seg[config['files']['segmentation_file']['column_names']['major_cn']],
        minor_cn=seg[config['files']['segmentation_file']['column_names']['minor_cn']],
        cell_frac=seg[config['files']['segmentation_file']['column_names']['cell_frac']],
        gene=seg[config['files']['segmentation_file']['column_names']['gene']]))
# if only integer copy number information is present
elif (options.seg_file and
        (config['files']['segmentation_file']['column_names']['major_cn'] in seg.columns and
        config['files']['segmentation_file']['column_names']['minor_cn'] in seg.columns)):
    source_seg = ColumnDataSource(data=dict(
        chrom=seg[config['files']['segmentation_file']['column_names']['chromosome']],
        start=seg[config['files']['segmentation_file']['column_names']['start']],
        end=seg[config['files']['segmentation_file']['column_names']['end']],
        logFC=seg[config['files']['segmentation_file']['column_names']['log2FC']],
        prob_start=seg.prob_start,
        prob_end=seg.prob_end,
        genome_cumsum_start=seg.genome_cumsum_start,
        genome_cumsum_end=seg.genome_cumsum_end,
        major_cn=seg[config['files']['segmentation_file']['column_names']['major_cn']],
        minor_cn=seg[config['files']['segmentation_file']['column_names']['minor_cn']],
        gene=seg[config['files']['segmentation_file']['column_names']['gene']]))
# if only clonality information is present
elif (options.seg_file and
    config['files']['segmentation_file']['column_names']['cell_frac'] in seg.columns):
    source_seg = ColumnDataSource(data=dict(
        chrom=seg[config['files']['segmentation_file']['column_names']['chromosome']],
        start=seg[config['files']['segmentation_file']['column_names']['start']],
        end=seg[config['files']['segmentation_file']['column_names']['end']],
        logFC=seg[config['files']['segmentation_file']['column_names']['log2FC']],
        prob_start=seg.prob_start,
        prob_end=seg.prob_end,
        genome_cumsum_start=seg.genome_cumsum_start,
        genome_cumsum_end=seg.genome_cumsum_end,
        cell_frac=round(seg[config['files']['segmentation_file']['column_names']['cell_frac']],2),
        gene=seg[config['files']['segmentation_file']['column_names']['gene']]))
elif (options.seg_file):
    source_seg = ColumnDataSource(data=dict(
        chrom=seg[config['files']['segmentation_file']['column_names']['chromosome']],
        start=seg[config['files']['segmentation_file']['column_names']['start']],
        end=seg[config['files']['segmentation_file']['column_names']['end']],
        logFC=seg[config['files']['segmentation_file']['column_names']['log2FC']],
        prob_start=seg.prob_start,
        prob_end=seg.prob_end,
        genome_cumsum_start=seg.genome_cumsum_start,
        genome_cumsum_end=seg.genome_cumsum_end,
        gene=seg[config['files']['segmentation_file']['column_names']['gene']]))


if (options.gene_file):
    source_gene = ColumnDataSource(data=dict(
        chrom=gen[config['files']['gene_file']['column_names']['chromosome']],
        start=gen[config['files']['gene_file']['column_names']['start']],
        end=gen[config['files']['gene_file']['column_names']['end']],
        logFC=gen[config['files']['gene_file']['column_names']['log2FC']],
        prob_start=gen.prob_start,
        prob_end=gen.prob_end,
        genome_cumsum_start=gen.genome_cumsum_start,
        genome_cumsum_end=gen.genome_cumsum_end,
        gene=gen[config['files']['gene_file']['column_names']['gene']]))

if (options.seg_blacklist and config['plots']['logFC_genome_plot']['artifact_mask']['visibility'] == "on"):
    source_seg_blacklist = ColumnDataSource(data=dict(
        chrom=seg_blacklist['chromosome'],
        start=seg_blacklist['start'],
        end=seg_blacklist['end'],
        genome_cumsum_start=seg_blacklist['genome_cumsum_start'],
        genome_cumsum_end=seg_blacklist['genome_cumsum_end']))

if (options.annot_file):
    source_gene_track = ColumnDataSource(data=dict(
        chrom=gene_track[config['files']['annotation_file']['column_names']['chromosome']],
        start=gene_track[config['files']['annotation_file']['column_names']['exon_start']],
        end=gene_track[config['files']['annotation_file']['column_names']['exon_end']],
        txid=gene_track[config['files']['annotation_file']['column_names']['tx_id']],
        exonNum=gene_track[config['files']['annotation_file']['column_names']['exon_number']],
        genome_cumsum_start=gene_track.genome_cumsum_start,
        genome_cumsum_end=gene_track.genome_cumsum_end,
        gene=gene_track[config['files']['annotation_file']['column_names']['gene']]))

if (options.vcf_file and not df_vaf.empty):
    source_vaf = ColumnDataSource(data=dict(
        chrom=df_vaf.chromosome,
        start=df_vaf.start,
        genome_cumsum=df_vaf.genome_cumsum,
        AF=df_vaf.AF,
        DP=df_vaf.DP,
        Ref=df_vaf.ref.astype(str),
        Alt=df_vaf.alt.astype(str)))

    if (options.seg_file):
        if (config['files']['segmentation_file']['column_names']['start'] == "start"):
            source_vaf_seg = ColumnDataSource(data=dict(
                chrom=df_vaf_seg[config['files']['segmentation_file']['column_names']['chromosome']],
                start=df_vaf_seg.start_y,
                end=df_vaf_seg[config['files']['segmentation_file']['column_names']['end']],
                AF=df_vaf_seg.AF,
                genome_cumsum_start=df_vaf_seg.genome_cumsum_start,
                genome_cumsum_end=df_vaf_seg.genome_cumsum_end))
        else:
            source_vaf_seg = ColumnDataSource(data=dict(
                chrom=df_vaf_seg[config['files']['segmentation_file']['column_names']['chromosome']],
                start=df_vaf_seg[config['files']['segmentation_file']['column_names']['start']],
                end=df_vaf_seg[config['files']['segmentation_file']['column_names']['end']],
                AF=df_vaf_seg.AF,
                genome_cumsum_start=df_vaf_seg.genome_cumsum_start,
                genome_cumsum_end=df_vaf_seg.genome_cumsum_end))

    if (options.gene_file):
        source_vaf_gene = ColumnDataSource(data=dict(
            chrom=df_vaf_gene.chromosome,
            start=df_vaf_gene.gene_start,
            end=df_vaf_gene.gene_end,
            AF=df_vaf_gene.AF,
            genome_cumsum_start=df_vaf_gene.genome_cumsum_start,
            genome_cumsum_end=df_vaf_gene.genome_cumsum_end))

TOOLTIPS_GENE_TRACK = [
    #("Index", "$index"),
    ("Chr", "@chrom"),
    ("Start - End", "@start - @end"),
    ("Gene", "@gene"),
    ("TxID | ExonNum", "@txid | @exonNum")
    #("Start", "@start"),
    #("End", "@end"),
    #("ExonNum", "@exonNum")
]

TOOLTIPS_INT_CN = [
    #("Index", "$index"),
    ("Cell Fraction", "@cell_frac{0.00}"),
    # ("Start - End", "@start - @end"),
    # ("Gene", "@gene"),
    # ("TxID | ExonNum", "@txid | @exonNum")
    #("Start", "@start"),
    #("End", "@end"),
    #("ExonNum", "@exonNum")
]

TOOLTIPS = [
    #("Index", "$index"),
    ("Chr", "@chrom"),
    ("Start - End", "@start - @end"),
    ("Gene", "@gene"),
    #("Start", "@start"),
    #("End", "@end"),
    ("Log2(FC)", "@logFC")
]

if (options.vcf_file and not df_vaf.empty):
    TOOLTIPS_VAF = [
        ("Chr", "@chrom"),
        ("Pos", "@start"),
        ("Ref | Alt", "@Ref | @Alt"),
        #("Alt", "@Alt"),
        ("DP", "@DP"),
        ("VAF", "@AF")
    ]

logFC = figure(plot_width=config['plots']['logFC_ind_plot']['width'],
               plot_height=config['plots']['logFC_ind_plot']['height'],
               tooltips=TOOLTIPS,
               tools=config['plots']['plot_tools'],
               output_backend=config['plots']['logFC_ind_plot']['output_backend'],
               active_scroll=config['plots']['logFC_ind_plot']['active_scroll'],
               active_tap="auto",
               title=config['plots']['logFC_ind_plot']['title'],
               #x_range=DataRange1d(bounds='auto'),
               x_range=DataRange1d(bounds=(min(data['ind']) - 0.05 * max(data['ind']),
                                           max(data['ind']) + 0.05 * max(data['ind']))),
               y_range=DataRange1d(bounds=(min(data[config['files']['ratio_file']['column_names']['log2FC']]) - 1,
                                           max(data[config['files']['ratio_file']['column_names']['log2FC']]) + 1)))

if (options.gene_file and config['plots']['logFC_ind_plot']['gene_markers']['visibility'] == "on"):
    logFC.quad(top="logFC",
               bottom=0,
               left="prob_start",
               right="prob_end",
               color=config['plots']['logFC_ind_plot']['gene_markers']['color'],
               source=source_gene)

logFC.circle("ind", "logFC",
             source=source,
             size="weight",
             line_color=config['plots']['logFC_ind_plot']['point_line_color'],
             fill_color=factor_cmap("label", palette=[config['plots']['logFC_ind_plot']['point_on_target_color'],
                                                      config['plots']['logFC_ind_plot']['point_off_target_color']],
                                    factors=["Target", "Antitarget"]))

if (options.seg_file and config['plots']['logFC_ind_plot']['segment_markers']['visibility'] == "on"):
    logFC.segment(x0="prob_start",
                  y0="logFC",
                  x1="prob_end",
                  y1="logFC",
                  color=config['plots']['logFC_ind_plot']['segment_markers']['segment_line_color'],
                  line_width=config['plots']['logFC_ind_plot']['segment_markers']['segment_line_width'],
                  level="annotation",
                  alpha=config['plots']['logFC_ind_plot']['segment_markers']['segment_line_alpha'],
                  source=source_seg)  # segmentation line

logFC.xaxis.axis_label = config['plots']['logFC_ind_plot']['x_axis_label']
logFC.yaxis.axis_label = config['plots']['logFC_ind_plot']['y_axis_label']
logFC.xaxis.visible = config['plots']['logFC_ind_plot']['x_axis_label_visibility'] == "on"
logFC.yaxis.visible = config['plots']['logFC_ind_plot']['y_axis_label_visibility'] == "on"
logFC.xgrid.grid_line_color = None

#####################

logFC_genome = figure(plot_width=config['plots']['logFC_genome_plot']['width'],
                      plot_height=config['plots']['logFC_genome_plot']['height'],
                      tooltips=TOOLTIPS,
                      tools=config['plots']['plot_tools'],
                      output_backend=config['plots']['logFC_genome_plot']['output_backend'],
                      active_scroll=config['plots']['logFC_genome_plot']['active_scroll'],
                      active_tap="auto",
                      x_range=DataRange1d(bounds=(min(data.genome_cumsum) - 0.05 * max(data.genome_cumsum),
                                                  max(data.genome_cumsum) + 0.05 * max(data.genome_cumsum))),
                      y_range=DataRange1d(
                          bounds=(min(data[config['files']['ratio_file']['column_names']['log2FC']]) - 1,
                                  max(data[config['files']['ratio_file']['column_names']['log2FC']]) + 1)),
                      title=config['plots']['logFC_genome_plot']['title'] + " (" + title_name + ")" + "\t\tGender: " +
                            str(options.gender) + "\t\tPurity: " + str(options.purity) +
                            "\t\tPloidy: " + str(options.ploidy))

if (options.gene_file and config['plots']['logFC_genome_plot']['gene_markers']['visibility'] == "on"):
    logFC_genome.quad(top="logFC",
                      bottom=0,
                      left="genome_cumsum_start",
                      right="genome_cumsum_end",
                      color=config['plots']['logFC_genome_plot']['gene_markers']['color'],
                      source=source_gene)

logFC_genome.circle("genome_cumsum", "logFC",
                    source=source,
                    size="weight",
                    line_color=config['plots']['logFC_genome_plot']['point_line_color'],
                    fill_color=factor_cmap("label",
                                           palette=[config['plots']['logFC_genome_plot']['point_on_target_color'],
                                                    config['plots']['logFC_genome_plot']['point_off_target_color']],
                                           factors=["Target", "Antitarget"]))

if (options.seg_file and config['plots']['logFC_genome_plot']['segment_markers']['visibility'] == "on"):
    logFC_genome.segment(x0="genome_cumsum_start",
                         y0="logFC",
                         x1="genome_cumsum_end",
                         y1="logFC",
                         color=config['plots']['logFC_genome_plot']['segment_markers']['segment_line_color'],
                         line_width=config['plots']['logFC_genome_plot']['segment_markers']['segment_line_width'],
                         level="annotation",
                         alpha=config['plots']['logFC_genome_plot']['segment_markers']['segment_line_alpha'],
                         source=source_seg)  # segmentation line

# print(options.bed_blacklist)
if (options.seg_blacklist and config['plots']['logFC_genome_plot']['artifact_mask']['visibility'] == "on"):
    # print("Entered!")
    logFC_genome.quad(top=max(data[config['files']['ratio_file']['column_names']['log2FC']]),
                      bottom=min(data[config['files']['ratio_file']['column_names']['log2FC']]),
                      left="genome_cumsum_start",
                      right="genome_cumsum_end",
                      color=config['plots']['logFC_genome_plot']['artifact_mask']['color'],
                      alpha=config['plots']['logFC_genome_plot']['artifact_mask']['alpha'],
                      level="underlay",
                      source=source_seg_blacklist)


logFC_genome.xaxis.axis_label = config['plots']['logFC_genome_plot']['x_axis_label']
logFC_genome.yaxis.axis_label = config['plots']['logFC_genome_plot']['y_axis_label']
logFC_genome.xaxis.visible = config['plots']['logFC_genome_plot']['x_axis_label_visibility'] == "on"
logFC_genome.yaxis.visible = config['plots']['logFC_genome_plot']['y_axis_label_visibility'] == "on"
logFC_genome.xgrid.grid_line_color = None

#############
int_cn_flag = 0

if(options.seg_file):
    int_cn_visibitlity = config['plots']['int_cn_plot']['visibility'] == "on" and \
            config['files']['segmentation_file']['column_names']['major_cn'] in seg.columns and \
            config['files']['segmentation_file']['column_names']['minor_cn'] in seg.columns
    # print(int_cn_visibitlity)
    # print(seg.columns)

    int_cn_flag = options.seg_file and \
            config['plots']['int_cn_plot']['visibility'] == "on" and \
            ((config['files']['segmentation_file']['column_names']['major_cn'] in seg.columns and
            config['files']['segmentation_file']['column_names']['minor_cn'] in seg.columns) or
            config['files']['segmentation_file']['column_names']['cell_frac'] in seg.columns)


if (options.seg_file and int_cn_flag):
    int_cn_genome = figure(plot_width=config['plots']['int_cn_plot']['width'],
                      plot_height=config['plots']['int_cn_plot']['height'],
                      tooltips=TOOLTIPS_INT_CN,
                      tools=config['plots']['plot_tools'],
                      output_backend=config['plots']['int_cn_plot']['output_backend'],
                      active_scroll=config['plots']['int_cn_plot']['active_scroll'],
                      active_tap="auto",
                      x_range=logFC_genome.x_range,
                      # y_range=DataRange1d(
                      #     bounds=(min(seg[config['files']['segmentation_file']['column_names']['minor_cn']]) - 1,
                      #             max(seg[config['files']['segmentation_file']['column_names']['major_cn']]) + 1)),
                      title=config['plots']['int_cn_plot']['title'])

    if (config['files']['segmentation_file']['column_names']['major_cn'] in seg.columns and
      config['files']['segmentation_file']['column_names']['minor_cn'] in seg.columns):

        int_cn_genome.segment(x0="genome_cumsum_start",
                             y0="minor_cn",
                             x1="genome_cumsum_end",
                             y1="minor_cn",
                             color=config['plots']['int_cn_plot']['segment_markers']['minor_segment_line_color'],
                             line_width=config['plots']['int_cn_plot']['segment_markers']['minor_segment_line_width'],
                             level="annotation",
                             alpha=config['plots']['int_cn_plot']['segment_markers']['minor_segment_line_alpha'],
                             source=source_seg)  # segmentation line

        int_cn_genome.segment(x0="genome_cumsum_start",
                             y0="major_cn",
                             x1="genome_cumsum_end",
                             y1="major_cn",
                             color=config['plots']['int_cn_plot']['segment_markers']['major_segment_line_color'],
                             line_width=config['plots']['int_cn_plot']['segment_markers']['major_segment_line_width'],
                             level="annotation",
                             alpha=config['plots']['int_cn_plot']['segment_markers']['major_segment_line_alpha'],
                             source=source_seg)  # segmentation line

    if (config['files']['segmentation_file']['column_names']['cell_frac'] in seg.columns):
        int_cn_genome.quad(top=max(seg[config['files']['segmentation_file']['column_names']['major_cn']]),
               bottom=min(seg[config['files']['segmentation_file']['column_names']['minor_cn']]),
               left="genome_cumsum_start",
               right="genome_cumsum_end",
               color=config['plots']['int_cn_plot']['cell_frac_color'],
               alpha="cell_frac",
               source=source_seg)


    int_cn_genome.xaxis.axis_label = config['plots']['int_cn_plot']['x_axis_label']
    int_cn_genome.yaxis.axis_label = config['plots']['int_cn_plot']['y_axis_label']
    int_cn_genome.xaxis.visible = config['plots']['int_cn_plot']['x_axis_label_visibility'] == "on"
    int_cn_genome.yaxis.visible = config['plots']['int_cn_plot']['y_axis_label_visibility'] == "on"
    int_cn_genome.xgrid.grid_line_color = None


#############
if (options.annot_file):
    logFC_genome_gene_track = figure(plot_width=config['plots']['annotation_plot']['width'],
                                     plot_height=config['plots']['annotation_plot']['height'],
                                     x_range=logFC_genome.x_range,
                                     tooltips=TOOLTIPS_GENE_TRACK,
                                     tools=config['plots']['plot_tools'],
                                     output_backend=config['plots']['annotation_plot']['output_backend'],
                                     active_scroll=config['plots']['annotation_plot']['active_scroll'],
                                     active_tap="auto",
                                     y_range=DataRange1d(bounds=(-2, 2)),
                                     title=config['plots']['annotation_plot']['title'])

    logFC_genome_gene_track.segment(x0="genome_cumsum_start",
                                    y0=0,
                                    x1="genome_cumsum_end",
                                    y1=0,
                                    color=config['plots']['annotation_plot']['line_color'],
                                    line_width=config['plots']['annotation_plot']['line_width'],
                                    level="annotation",
                                    alpha=config['plots']['annotation_plot']['line_alpha'],
                                    source=source_gene_track)

    logFC_genome_gene_track.yaxis.visible = False
    logFC_genome_gene_track.xaxis.visible = False
    logFC_genome_gene_track.xgrid.grid_line_color = None
    logFC_genome_gene_track.ygrid.grid_line_color = None

#############
if (options.vcf_file and not df_vaf.empty):
    VAF_genome = figure(plot_width=config['plots']['vaf_plot']['width'],
                        plot_height=config['plots']['vaf_plot']['height'],
                        x_range=logFC_genome.x_range,  # y_range = -0.7,
                        tooltips=TOOLTIPS_VAF,
                        tools=config['plots']['plot_tools'],
                        output_backend=config['plots']['vaf_plot']['output_backend'],
                        active_scroll=config['plots']['vaf_plot']['active_scroll'],
                        active_tap="auto",
                        y_range=DataRange1d(bounds=(-0.05, 1.05)),
                        title=config['plots']['vaf_plot']['title'])

    VAF_genome.circle("genome_cumsum", "AF",
                      source=source_vaf,
                      size=config['plots']['vaf_plot']['point_size'],
                      line_color=config['plots']['vaf_plot']['point_line_color'],
                      fill_color=config['plots']['vaf_plot']['point_color'])

    if (options.seg_file and config['plots']['vaf_plot']['segment_markers']['visibility'] == "on"):
        VAF_genome.segment(x0="genome_cumsum_start",
                           y0="AF",
                           x1="genome_cumsum_end",
                           y1="AF",
                           color=config['plots']['vaf_plot']['segment_markers']['segment_line_color'],
                           line_width=config['plots']['vaf_plot']['segment_markers']['segment_line_width'],
                           level="annotation",
                           alpha=config['plots']['vaf_plot']['segment_markers']['segment_line_alpha'],
                           source=source_vaf_seg)  # segmentation line

    if (options.gene_file and config['plots']['vaf_plot']['gene_markers']['visibility'] == "on"):
        VAF_genome.segment(x0="genome_cumsum_start",
                           y0="AF",
                           x1="genome_cumsum_end",
                           y1="AF",
                           color=config['plots']['vaf_plot']['gene_markers']['gene_line_color'],
                           line_width=config['plots']['vaf_plot']['gene_markers']['gene_line_width'],
                           level="annotation",
                           alpha=config['plots']['vaf_plot']['gene_markers']['gene_line_alpha'],
                           source=source_vaf_gene)  # gene average line

    VAF_genome.xaxis.axis_label = config['plots']['vaf_plot']['x_axis_label']
    VAF_genome.yaxis.axis_label = config['plots']['vaf_plot']['y_axis_label']
    VAF_genome.xaxis.visible = config['plots']['vaf_plot']['x_axis_label_visibility'] == "on"
    VAF_genome.yaxis.visible = config['plots']['vaf_plot']['y_axis_label_visibility'] == "on"
    VAF_genome.xgrid.grid_line_color = None

    logFC_genome.min_border_bottom = 0
    VAF_genome.min_border_top = 0

# Setting up taptool for the link out to UCSC Genome Browser
url = "https://genome.ucsc.edu/cgi-bin/hgTracks?db=" + config['files']['genome_build'] + "&position=@chrom:@start-@end"
taptool = logFC.select(type=TapTool)
taptool.callback = OpenURL(url=url)
taptool = logFC_genome.select(type=TapTool)
taptool.callback = OpenURL(url=url)

chr_boundary = pd.concat([data[config['files']['ratio_file']['column_names']['chromosome']],
                          data.ind],
                         axis=1).groupby([config['files']['ratio_file']['column_names']['chromosome']],
                                         as_index=False).max()
# print(chr_boundary)
draw_chr_boundary(logFC, chr_boundary, genome=False, vaf=False)

chr_boundary = pd.concat([data[config['files']['ratio_file']['column_names']['chromosome']],
                          data.genome_cumsum],
                         axis=1).groupby([config['files']['ratio_file']['column_names']['chromosome']],
                                         as_index=False).max()
draw_chr_boundary(logFC_genome, chr_boundary, genome=True, vaf=False)

if (options.vcf_file and not df_vaf.empty):
    chr_boundary = pd.concat([data[config['files']['ratio_file']['column_names']['chromosome']],
                              data.genome_cumsum],
                             axis=1).groupby([config['files']['ratio_file']['column_names']['chromosome']],
                                             as_index=False).max()
    draw_chr_boundary(VAF_genome, chr_boundary, genome=True, vaf=True)

    taptool = VAF_genome.select(type=TapTool)
    taptool.callback = OpenURL(url=url)

# Setting up data table
columns = [
    TableColumn(field="chrom", title="Chr"),
    TableColumn(field="start", title="Start"),
    TableColumn(field="end", title="End"),
    TableColumn(field="gene", title="Gene"),
    TableColumn(field="logFC", title="Log2(FC)", formatter=NumberFormatter(format="0.00")),
]

data_table_cnr = DataTable(source=source,
                           columns=columns,
                           fit_columns=True,
                           width=config['plots']['ratio_table']['width'],
                           height=config['plots']['ratio_table']['height'])

if (options.gene_file):
    data_table_gene = DataTable(source=source_gene,
                                columns=columns,
                                fit_columns=True,
                                width=config['plots']['gene_table']['width'],
                                height=config['plots']['gene_table']['height'])

    selection_options = ["---", "All", "All Amplification and Loss", "Amplification", "Gain", "Loss", "Deep Loss",
                         "---", config['files']['ratio_file']['off_target_label']]
    selection_options.extend(np.unique(gen[config['files']['gene_file']['column_names']['gene']].tolist()))

    select = MultiSelect(title="Gene filter", options=selection_options, max_width=300, size=15)

    filteredSource = ColumnDataSource(data=dict(chrom=[], start=[], end=[], gene=[], logFC=[]))

    data_table = DataTable(source=filteredSource,
                           columns=columns,
                           fit_columns=True,
                           width=config['plots']['filtered_table']['width'],
                           height=config['plots']['filtered_table']['height'])

    callback_filter = CustomJS(args=dict(source=source,
                                         source_copy=source_copy,
                                         source_gene=source_gene,
                                         filteredSource=filteredSource,
                                         data_table=data_table_gene,
                                         loss_threshold=config['files']['gene_file']['loss_threshold'],
                                         deep_loss_threshold = config['files']['gene_file']['deep_loss_threshold'],
                                         gain_threshold=config['files']['gene_file']['gain_threshold'],
                                         amp_threshold=config['files']['gene_file']['amp_threshold']), code="""
    var data_gene = source_gene.data;
    var f = cb_obj.value;
    var d2 = filteredSource.data;
    d2['chrom']=[];
    d2['start']=[];
    d2['end']=[];
    d2['gene']=[];
    d2['logFC']=[];
    var thresh = 0;

    for(i = 0; i < data_gene['chrom'].length;i++){

    if(f == "Deep Loss"){
        thresh = deep_loss_threshold
        if(data_gene['logFC'][i] <= thresh){
            d2['chrom'].push(data_gene['chrom'][i])
            d2['start'].push(data_gene['start'][i])
            d2['end'].push(data_gene['end'][i])
            d2['gene'].push(data_gene['gene'][i])
            d2['logFC'].push(data_gene['logFC'][i])
        }
    }else if(f == "Loss"){
        thresh = loss_threshold
        if(data_gene['logFC'][i] <= thresh & data_gene['logFC'][i] > deep_loss_threshold){
            d2['chrom'].push(data_gene['chrom'][i])
            d2['start'].push(data_gene['start'][i])
            d2['end'].push(data_gene['end'][i])
            d2['gene'].push(data_gene['gene'][i])
            d2['logFC'].push(data_gene['logFC'][i])
        }
    }else if(f == "Amplification"){
        thresh = amp_threshold
        if(data_gene['logFC'][i] >= thresh){
            d2['chrom'].push(data_gene['chrom'][i])
            d2['start'].push(data_gene['start'][i])
            d2['end'].push(data_gene['end'][i])
            d2['gene'].push(data_gene['gene'][i])
            d2['logFC'].push(data_gene['logFC'][i])
        }
    }else if(f == "Gain"){
        if(data_gene['logFC'][i] >= gain_threshold & data_gene['logFC'][i] < amp_threshold){
            d2['chrom'].push(data_gene['chrom'][i])
            d2['start'].push(data_gene['start'][i])
            d2['end'].push(data_gene['end'][i])
            d2['gene'].push(data_gene['gene'][i])
            d2['logFC'].push(data_gene['logFC'][i])
        }
    }else if (f == "All"){
        thresh = 0
        d2['chrom'].push(data_gene['chrom'][i])
        d2['start'].push(data_gene['start'][i])
        d2['end'].push(data_gene['end'][i])
        d2['gene'].push(data_gene['gene'][i])
        d2['logFC'].push(data_gene['logFC'][i])

    }else if (f == "All Amplification and Loss"){
        thresh_amp = amp_threshold
        thresh_del = loss_threshold
        if(data_gene['logFC'][i] >= thresh_amp | data_gene['logFC'][i] <= thresh_del){
            d2['chrom'].push(data_gene['chrom'][i])
            d2['start'].push(data_gene['start'][i])
            d2['end'].push(data_gene['end'][i])
            d2['gene'].push(data_gene['gene'][i])
            d2['logFC'].push(data_gene['logFC'][i])
        }

    }else if (f.includes(data_gene['gene'][i])){
        d2['chrom'].push(data_gene['chrom'][i])
        d2['start'].push(data_gene['start'][i])
        d2['end'].push(data_gene['end'][i])
        d2['gene'].push(data_gene['gene'][i])
        d2['logFC'].push(data_gene['logFC'][i])

    }
    }

    var data = source_copy.data;
    var d3 = {};
    d3['chrom']=[];
    d3['start']=[];
    d3['genome_cumsum']=[];
    d3['end']=[];
    d3['gene']=[];
    d3['logFC']=[];
    d3['ind']=[];
    d3['label']=[];
    d3['weight']=[];


    for(i = 0; i < data['chrom'].length;i++){
    if (f == "All" ){
        thresh = 0
        d3['chrom'].push(data['chrom'][i])
        d3['start'].push(data['start'][i])
        d3['end'].push(data['end'][i])
        d3['gene'].push(data['gene'][i])
        d3['logFC'].push(data['logFC'][i])
        d3['genome_cumsum'].push(data['genome_cumsum'][i])
        d3['ind'].push(data['ind'][i])
        d3['label'].push(data['label'][i])
        d3['weight'].push(data['weight'][i])

    }else if(f == "Loss" | f == "Amplification" | f == "All Amplification and Loss" | f == "Gain" | f == "Deep Loss"){
        if(d2['gene'].includes(data['gene'][i])){
            d3['chrom'].push(data['chrom'][i])
            d3['start'].push(data['start'][i])
            d3['end'].push(data['end'][i])
            d3['gene'].push(data['gene'][i])
            d3['logFC'].push(data['logFC'][i])
            d3['genome_cumsum'].push(data['genome_cumsum'][i])
            d3['ind'].push(data['ind'][i])
            d3['label'].push(data['label'][i])
            d3['weight'].push(data['weight'][i])
        }
    } else if (f.includes(data['gene'][i])){
        d3['chrom'].push(data['chrom'][i])
        d3['start'].push(data['start'][i])
        d3['end'].push(data['end'][i])
        d3['gene'].push(data['gene'][i])
        d3['logFC'].push(data['logFC'][i])
        d3['genome_cumsum'].push(data['genome_cumsum'][i])
        d3['ind'].push(data['ind'][i])
        d3['label'].push(data['label'][i])
        d3['weight'].push(data['weight'][i])

    } 
    }
    
    if (!d3['chrom'].length){
        d3 = data
    }
    
    source.data = d3
    filteredSource.change.emit()
    // trigger change on datatable
    data_table.change.emit()

    """)

    select.js_on_change('value', callback_filter)

    table_data_div = Div(text="<b>Filtered Data</b>")
    table_gene_div = Div(text="<b>Gene Level Calling</b>")

table_bin_div = Div(text="<b>Bin Data</b>")

if (options.vcf_file and options.annot_file and options.gene_file and not df_vaf.empty):
    if (int_cn_flag):
        plots = gridplot([[logFC_genome], [logFC_genome_gene_track], [VAF_genome], [int_cn_genome],
                          [logFC]], toolbar_location='left', merge_tools=True)
    else :
        plots = gridplot([[logFC_genome], [logFC_genome_gene_track], [VAF_genome], [logFC]], toolbar_location='left',
                         merge_tools=True)
    fig_datatables = layout(
        row(column(table_gene_div, data_table_gene), column(table_data_div, data_table), column(select)))
    final_fig = layout(children=[[plots], [fig_datatables]])

elif (options.vcf_file and options.annot_file and not df_vaf.empty):
    if (int_cn_flag):
        plots = gridplot([[logFC_genome], [logFC_genome_gene_track], [VAF_genome], [int_cn_genome],
                          [logFC]], toolbar_location='left', merge_tools=True)
    else:
        plots = gridplot([[logFC_genome], [logFC_genome_gene_track], [VAF_genome], [logFC]], toolbar_location='left',
                         merge_tools=True)
    fig_datatables = layout(
        row(column(table_bin_div, data_table_cnr)))
    final_fig = layout(children=[[plots], [fig_datatables]])

elif (options.vcf_file and options.gene_file and not df_vaf.empty):
    if (int_cn_flag):
        plots = gridplot([[logFC_genome], [VAF_genome], [int_cn_genome], [logFC]], toolbar_location='left',
                     merge_tools=True)
    else:
        plots = gridplot([[logFC_genome], [VAF_genome], [logFC]], toolbar_location='left',
                         merge_tools=True)
    fig_datatables = layout(
        row(column(table_gene_div, data_table_gene), column(table_data_div, data_table), column(select)))
    final_fig = layout(children=[[plots], [fig_datatables]])

elif (options.annot_file and options.gene_file):
    if (int_cn_flag):
        plots = gridplot([[logFC_genome], [logFC_genome_gene_track], [int_cn_genome], [logFC]],
                         toolbar_location='left', merge_tools=True)
    else:
        plots = gridplot([[logFC_genome], [logFC_genome_gene_track], [logFC]], toolbar_location='left',
                         merge_tools=True)
    fig_datatables = layout(
        row(column(table_gene_div, data_table_gene), column(table_data_div, data_table), column(select)))
    final_fig = layout(children=[[plots], [fig_datatables]])

elif (options.gene_file):
    if (int_cn_flag):
        plots = gridplot([[logFC_genome], [int_cn_genome], [logFC]], toolbar_location='left',
                         merge_tools=True)
    else:
        plots = gridplot([[logFC_genome], [logFC]], toolbar_location='left',
                         merge_tools=True)
    fig_datatables = layout(
        row(column(table_gene_div, data_table_gene), column(table_data_div, data_table), column(select)))
    final_fig = layout(children=[[plots], [fig_datatables]])

elif (options.annot_file):
    if (int_cn_flag):
        plots = gridplot([[logFC_genome], [logFC_genome_gene_track], [int_cn_genome], [logFC]],
                         toolbar_location='left', merge_tools=True)
    else:
        plots = gridplot([[logFC_genome], [logFC_genome_gene_track], [logFC]], toolbar_location='left',
                         merge_tools=True)
    fig_datatables = layout(
        row(column(table_bin_div, data_table_cnr)))
    final_fig = layout(children=[[plots], [fig_datatables]])

elif (options.vcf_file and not df_vaf.empty):
    if (int_cn_flag):
        plots = gridplot([[logFC_genome], [VAF_genome], [int_cn_genome], [logFC]], toolbar_location='left',
                     merge_tools=True)
    else:
        plots = gridplot([[logFC_genome], [VAF_genome], [logFC]], toolbar_location='left',
                         merge_tools=True)
    fig_datatables = layout(
        row(column(table_bin_div, data_table_cnr)))
    final_fig = layout(children=[[plots], [fig_datatables]])

else:
    if (options.seg_file and int_cn_flag):
        plots = gridplot([[logFC_genome], [int_cn_genome], [logFC]], toolbar_location='left',
                         merge_tools=True)
    else:
        plots = gridplot([[logFC_genome], [logFC]], toolbar_location='left',
                         merge_tools=True)
    fig_datatables = layout(
        row(column(table_bin_div, data_table_cnr)))
    final_fig = layout(children=[[plots], [fig_datatables]])

if (config['plots']['bokeh_js_css_code'] == "INLINE"):
    save(final_fig, resources=INLINE)
else:
    save(final_fig)
# show(final_fig)
logging.info("Done!")
