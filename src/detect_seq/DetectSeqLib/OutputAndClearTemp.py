# _*_ coding: UTF-8 _*_

import logging
import sys
import os

# Version information START ----------------------------------------------------
VERSION_INFO = \
"""

Author: MENG Howard

Version-01:
  2020-11-10 Merge files and remove temp files

E-Mail: meng_howard@126.com
"""
# Version information END ------------------------------------------------------


####################################################################################################
# Output part
####################################################################################################
def merge_split_files(split_file_dict,
                      key_order_list=None,
                      out_filename="stdout",
                      header_list=None,
                      in_sep="\t", out_sep="\t",
                      log_verbose=3,
                      return_col_index=None):
    """
    INPUT:
        <split_file_dict>
            dict, typical format like
            {
                'chr1': 'chr1.eqdA2zJkHCm0YW1o.AddBlockInfo',
                'chr19': 'chr19.ZIrU6xC2DFcBayE1.AddBlockInfo',
                'chr20': 'chr20.nK2goEB6xh9MzXpD.AddBlockInfo',
                'chr_name_order': ['chr1', 'chr19', 'chr20']
            }

        <key_order_list>
            list, if set 'None' as default, will merge all split files and make output


        <out_filename>
            str

        <header_list>
            list, length have to match with split files

        <sep>
            str, default '\t', can set by custom

    RETURN:
        0, everything is okay.

    INFO:
        Final date 2020-11-10
    """
    # log
    logging.basicConfig(level=(4 - log_verbose) * 10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w")

    # return col
    return_col_list = []

    # open output file
    if out_filename == "stdout":
        out_file = sys.stdout
    else:
        out_file = open(out_filename, "w")

    # init vars
    ignore_key_list = ["chr_name_order"]

    # check keys
    if key_order_list is None:
        key_order_list = split_file_dict.keys()

    # header output
    if header_list is not None:
        out_file.write(out_sep.join(map(str, header_list)) + "\n")

    # merge files
    logging.info("Starting to merge files...")

    for run_key in key_order_list:
        if run_key not in ignore_key_list:
            logging.debug("Merging files, processing on \n\t%s" % split_file_dict[run_key])

            with open(split_file_dict[run_key], "r") as run_file:
                for line in run_file:
                    if in_sep != out_sep:
                        line_list = line.strip().split(in_sep)
                        out_file.write(out_sep.join(line_list) + "\n")

                        if return_col_index is not None:
                            return_col_list.append(line_list[return_col_index])
                    else:
                        out_file.write(line)

                        if return_col_index is not None:
                            line_list = line.strip().split(in_sep)
                            return_col_list.append(line_list[return_col_index])

    # log
    logging.info("Merging files, done!")

    # close file
    if out_filename != "stdout":
        out_file.close()

    return 0, return_col_list


####################################################################################################
# clear temp
####################################################################################################
def clear_temp_files_by_dict(temp_file_dict, log_verbose=3):
    """
    INPUT:
        <temp_file_dict>
            dict, format like:
                {
                    'chr1': 'chr1.eqdA2zJkHCm0YW1o.AddBlockInfo',
                    'chr19': 'chr19.ZIrU6xC2DFcBayE1.AddBlockInfo',
                    'chr20': 'chr20.nK2goEB6xh9MzXpD.AddBlockInfo',
                    'chr_name_order': ['chr1', 'chr19', 'chr20']
                }

    RETURN:
        0, everything is okay.

    INFO:
        Final date 2020-11-10
    """

    # log
    logging.basicConfig(level=(4 - log_verbose) * 10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w")

    run_state = 0

    for run_key in temp_file_dict.keys():
        if type(temp_file_dict[run_key]) is str:
            if os.path.exists(temp_file_dict[run_key]):
                if os.path.isfile(temp_file_dict[run_key]):
                    try:
                        os.remove(temp_file_dict[run_key])
                    except:
                        logging.error("Removing error on \n\t%s" % temp_file_dict[run_key])
                        run_state = 1

    return run_state




















