#! /Users/meng/menghw_HD/miniconda3/bin/python
# _*_ coding: UTF-8 _*_

import argparse
import sys

import matplotlib
import numpy as np
import pandas as pd

from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

matplotlib.use('Agg')
np.set_printoptions(suppress=True)

# Version information START --------------------------------------------------
VERSION_INFO = """

Author: MENG Howard

Version-02:
    2022-02-13
        Support Python3
        Supoort for sgRNA alignment and TALE alignment

Version-01:
    2021-05-01
        Make plot for DdCBE TALE sequence

    2019-12-05 
        Plot Alignment Result Table (.art) files

E-Mail: meng_howard@126.com

"""
# Version information END ----------------------------------------------------


# Learning Part START --------------------------------------------------------
LEARNING_PART = \
    """

There will be a lot of parameters to set,

but default parameters can handle output figure in most cases.

"""

PIPELINE = \
    """

1. Load table 

2. sort table with key words

3. make plot

"""
# Learning Part END-----------------------------------------------------------

###############################################################################
# main part
###############################################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Plot (.art) file")

    parser.add_argument("-i", "--input_art",
                        help="Input alignment result table", required=True)

    parser.add_argument("--align_seq",
                        help="Alignment sequence", required=True)

    parser.add_argument("-k", "--sort_column_key",
                        help="Sort .art file by this columns order, use ',' as separator, default=None means don't "
                             "sort .art file",
                        default=None)

    parser.add_argument("-r", "--sort_ascend_state",
                        help="Ture means sort by ascending order, False means sort by descending order, "
                             "default=False,...",
                        default=None)

    parser.add_argument("-a", "--annotation_column_key",
                        help="The annotation info in output figure are default set as "
                             "align_total_mismatch,"
                             "align_degen_total_mismatch,"
                             "region_index,"
                             "align_coordinate,"
                             "align_strand",
                        default="align_total_mismatch,"
                                "align_degen_total_mismatch,"
                                "region_index,"
                                "align_coordinate,"
                                "align_strand")

    parser.add_argument("-o", "--out_figure",
                        help="Output figure filename", required=True)

    parser.add_argument("--out_figure_format",
                        help="Support 'pdf' and 'png' result, default=pdf", default="pdf")

    parser.add_argument("--out_figure_dpi",
                        help="Out figure dpi, default=96", default="96", type=int)

    ARGS = parser.parse_args()

    # ------------------------------------------------------------------------>>>>>>>
    # load sgRNA
    # ------------------------------------------------------------------------>>>>>>>
    sgRNA_full_seq = str(ARGS.align_seq).upper()
    sgRNA_model_seq = sgRNA_full_seq

    # ---------------------------------------------------------------->>>>>>>>
    # load art 
    # ---------------------------------------------------------------->>>>>>>>
    art_table = pd.read_csv(ARGS.input_art, sep="\t")

    # make coordinate
    align_coordinate_list = []
    for index in art_table.index:
        chr_name, align_start, align_end = art_table.loc[
            index, ["align_chr_name", "align_chr_start", "align_chr_end"]]
        align_coordinate_str = "{chr_name}:{chr_start}-{chr_end}".format(chr_name=chr_name,
                                                                         chr_start=align_start,
                                                                         chr_end=align_end)
        align_coordinate_list.append(align_coordinate_str)
    art_table["align_coordinate"] = align_coordinate_list

    # ---------------------------------------------------------------->>>>>>>>
    # sort art
    # ---------------------------------------------------------------->>>>>>>>
    if ARGS.sort_column_key is not None:
        sort_key_list = ARGS.sort_column_key.split(",")
        sort_key_check_state = all([x in art_table.columns for x in sort_key_list])

        if sort_key_check_state:
            if ARGS.sort_ascend_state is None:
                sort_ascending_list = [False] * len(sort_key_list)

            else:
                try:
                    sort_ascending_list = []
                    for sort_state_str in ARGS.sort_ascend_state.split(","):
                        if sort_state_str == "False":
                            sort_ascending_list.append(False)
                        else:
                            sort_ascending_list.append(True)

                except:
                    sys.stderr.write("WARNING:\n--sort_ascend_state format error, use default action!\n")
                    sort_ascending_list = [False] * len(sort_key_list)

                if len(sort_ascending_list) != len(sort_key_list):
                    sys.stderr.write("WARNING:\n--sort_ascend_state format error, use default action!\n")
                    sort_ascending_list = [False] * len(sort_key_list)

        art_table_filter = art_table.sort_values(by=sort_key_list, ascending=sort_ascending_list)

    else:
        art_table_filter = art_table

    # ---------------------------------------------------------------->>>>>>>>
    # annotation list
    # ---------------------------------------------------------------->>>>>>>>
    ann_colname_list = list(ARGS.annotation_column_key.split(","))

    ann_colname_check_state_list = []

    for ann_colname in ann_colname_list:
        if ann_colname in art_table.columns:
            ann_colname_check_state_list.append(True)
        else:
            ann_colname_check_state_list.append(False)

    ann_colname_check_state = all(ann_colname_check_state_list)

    if not ann_colname_check_state:
        sys.stderr.write("WARNING:\n--annotation_column_key format error, use default action!\n")
        ann_colname_list = ["align_total_mismatch",
                            "align_degen_total_mismatch",
                            "region_index",
                            "align_coordinate",
                            "align_strand"]

    # ---------------------------------------------------------------->>>>>>>>
    # plot region
    # ---------------------------------------------------------------->>>>>>>>

    # --------------------------------------------->>>>>
    # set all box position
    # --------------------------------------------->>>>>
    # set plot column num 
    box_column_num = len(sgRNA_model_seq) + len(ann_colname_list)
    box_row_num = len(art_table_filter) + 1

    panel_box_width = 0.3
    panel_box_height = 0.3
    panel_row_space = 0.03
    panel_col_space = 0.03

    panel_box_width_coef = [1] * len(sgRNA_model_seq)

    for colname in ann_colname_list:
        col_info_len_list = [len(str(col_info)) for col_info in art_table_filter[colname]]
        col_info_len_list.append(len(colname))
        col_width_coef = max(col_info_len_list) // 3 + 1
        panel_box_width_coef.append(col_width_coef)

    # set figure width and height 
    figure_width = panel_box_width * sum(panel_box_width_coef) + panel_col_space * (box_column_num - 1)
    figure_height = panel_box_height * box_row_num + panel_row_space * (box_row_num - 1)

    # calculate all box x vec
    box_x_left_vec = []
    box_x_right_vec = []

    x_left = x_right = 0
    for index, col_coef in enumerate(panel_box_width_coef):
        if index > 0:
            x_left += panel_box_width_coef[index - 1] * panel_box_width + panel_col_space
            x_right = x_left + panel_box_width_coef[index] * panel_box_width
        else:
            x_left = 0
            x_right = x_left + panel_box_width_coef[0] * panel_box_width

        box_x_left_vec.append(x_left)
        box_x_right_vec.append(x_right)

    # fix number 
    box_x_left_vec = np.array(box_x_left_vec)
    box_x_right_vec = np.array(box_x_right_vec)

    box_x_left_vec.round(decimals=5)
    box_x_right_vec.round(decimals=5)

    # --------------------------------------------->>>>>
    # plot alignment region
    # --------------------------------------------->>>>>
    box_border_plot_state = False

    # make box border 
    if box_border_plot_state:
        box_edgecolor = "#AAAAAA"
        box_linestyle = "-"
        box_linewidth = 2
    else:
        box_edgecolor = "#FFFFFF"
        box_linestyle = "None"
        box_linewidth = 0

    # base_color_dict = {"A": "#04E3E3", "T": "#F9B874", "C": "#B9E76B", "G": "#F53798", "N": "#AAAAAA", "-": "#AAAAAA"}
    base_color_dict = {"A": "#97C3E9", "T": "#F7B6B7", "C": "#A2DB69", "G": "#D868B9", "N": "#DBDBDB", "-": "#DBDBDB"}

    # set new figure 
    fig = plt.figure(figsize=(figure_width * 1.1, figure_height * 1.1))
    ax = fig.add_subplot(111, aspect="equal")
    plt.xlim([0, figure_width])
    plt.ylim([-figure_height, 0])
    plt.axis("off")

    # make plot 
    text_list = []
    patches = []

    current_y = 0
    # ---------------------------------------->>>>>
    # plot model and colname 
    # ---------------------------------------->>>>>
    # plot model
    col_index = 0
    box_y_btm = current_y - panel_box_height

    for index, rna_char in enumerate(sgRNA_model_seq):
        box_x_left = box_x_left_vec[index]
        box_width = box_x_right_vec[index] - box_x_left_vec[index]
        box_value = rna_char

        box_color = base_color_dict.get(box_value)
        if not box_color:
            box_fill = False
            box_color = "#FFFFFF"
        else:
            box_fill = True

        # plot sgRNA box
        patches.append(Rectangle(
            xy=(box_x_left, box_y_btm),
            width=panel_box_width_coef[index] * panel_box_width,
            height=panel_box_height,
            fill=box_fill,
            alpha=1,
            linestyle=box_linestyle,
            linewidth=box_linewidth,
            edgecolor=box_edgecolor,
            facecolor=box_color)
        )

        # add sgRNA text
        text_list.append((box_x_left + 0.5 * box_width, box_y_btm + 0.5 * panel_box_height, box_value, 16))

        # record index 
        col_index += 1

        # plot name
    for ann_colname in ann_colname_list:
        box_x_left = box_x_left_vec[col_index]
        box_width = box_x_right_vec[col_index] - box_x_left_vec[col_index]
        text_list.append((box_x_left + 0.5 * box_width, box_y_btm + 0.5 * panel_box_height, ann_colname, 10))
        col_index += 1

        # make cur_y
    current_y = current_y - panel_box_height - panel_row_space

    # ---------------------------------------->>>>>
    # plot remain info
    # ---------------------------------------->>>>>
    for index in art_table_filter.index:
        sgRNA_align, align_info, query_align = art_table_filter.loc[
            index, ["align_query_seq", "align_info_state", "align_target_seq"]]
        ann_value_list = list(art_table_filter.loc[index, ann_colname_list])

        col_index = 0
        sgRNA_gap_count = 0
        box_y_btm = current_y - panel_box_height

        for index, rna_char in enumerate(sgRNA_align):
            if rna_char != "-":
                sgRNA_model_char = sgRNA_model_seq[index - sgRNA_gap_count]
                if sgRNA_model_char == query_align[index]:
                    #                 box_value = r"$\bullet$"
                    box_value = r"$\cdot$"
                    box_fill = False
                    box_color = "#FFFFFF"
                else:
                    box_value = query_align[index]
                    box_fill = True
                    box_color = base_color_dict.get(box_value)

                    if not box_color:
                        box_color = "AAAAAA"
            else:
                sgRNA_gap_count += 1
                continue

            box_x_left = box_x_left_vec[index - sgRNA_gap_count]
            box_width = box_x_right_vec[index - sgRNA_gap_count] - box_x_left_vec[index - sgRNA_gap_count]

            # plot align seq box
            patches.append(Rectangle(
                xy=(box_x_left, box_y_btm),
                width=box_width,
                height=panel_box_height,
                fill=box_fill,
                alpha=1,
                linestyle=box_linestyle,
                linewidth=box_linewidth,
                edgecolor=box_edgecolor,
                facecolor=box_color)
            )

            # add sgRNA text
            text_list.append((box_x_left + 0.5 * box_width, box_y_btm + 0.5 * panel_box_height, str(box_value), 16))

            # record index 
            col_index += 1

            # plot name
        for ann_value in ann_value_list:
            box_x_left = box_x_left_vec[col_index]
            box_width = box_x_right_vec[col_index] - box_x_left_vec[col_index]
            try:
                ann_value_str = round(float(ann_value), 4)
            except:
                ann_value_str = str(ann_value)
            text_list.append((box_x_left + 0.5 * box_width, box_y_btm + 0.5 * panel_box_height, ann_value_str, 10))
            col_index += 1

        # make cur_y
        current_y = current_y - panel_box_height - panel_row_space

    # plot remain info
    ax.add_collection(PatchCollection(patches, match_original=True))

    # add text 
    for text_x, text_y, text_info, text_fontsize in text_list:
        plt.text(
            x=text_x,
            y=text_y,
            s=text_info,
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=text_fontsize,
            fontname="Arial"
        )

    fig.savefig(fname=ARGS.out_figure, bbox_inches='tight', dpi=ARGS.out_figure_dpi, format=ARGS.out_figure_format)

# 2022-02-13 by MENG Howard
