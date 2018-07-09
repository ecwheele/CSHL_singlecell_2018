#!/usr/bin/env python

"""
basic module to use as an example.
"""
import argparse
import pandas as pd
import collections
import tables
import scipy.io
import scipy.sparse
import csv
import numpy as np
from collections import defaultdict

def get_matrix_from_h5(filename, genome, master_barcode_file, barcode_mapper_output):
    """
    Returns matrix and attributes using cellranger's h5 accessor method.

    :param filename:
    :param genome:
    :return:
    """
    GeneBCMatrix = collections.namedtuple('GeneBCMatrix',
                                          ['gene_ids', 'gene_names',
                                           'barcodes', 'matrix'])

    with tables.open_file(filename, 'r') as f:
        try:
            group = f.get_node(f.root, genome)
        except tables.NoSuchNodeError:
            print("That genome does not exist in this file.")
            return None
        gene_ids = getattr(group, 'genes').read()
        gene_names = getattr(group, 'gene_names').read()
        barcodes = getattr(group, 'barcodes').read()
        data = getattr(group, 'data').read()
        # print('data', data)
        indices = getattr(group, 'indices').read()
        # print('indices', indices)
        indptr = getattr(group, 'indptr').read()
        # print('indptr', indptr)
        shape = getattr(group, 'shape').read()
        # print('shape', shape)


        matrix = scipy.sparse.csc_matrix((data, indices, indptr), shape=shape)

        if master_barcode_file is not None and barcode_mapper_output is not None:
            barcodes = map_barcodes(barcodes, master_barcode_file,
                                   barcode_mapper_output)
        elif master_barcode_file is not None:
            barcodes = unmap_barcodes(barcodes, master_barcode_file)


        return GeneBCMatrix(gene_ids, gene_names, barcodes, matrix)


def read_mtx_as_dataframe(mtx_file, columns_file, rows_file, master_barcode_file, barcode_mapper_output):
    """
    Reads a mtx file and returns a pandas dataframe
    :param mtx_file: sparse matrix
    :return df: Pandas.DataFrame()
    """
    mat = scipy.io.mmread(mtx_file)


    columns = [
        row[0] for row in csv.reader(open(columns_file), delimiter="\t")
    ]
    rows = [
        row[0] for row in csv.reader(open(rows_file), delimiter="\t")
    ]

    # if master_barcode_file is not None and barcode_mapper_output is not None:
    #     columns = map_barcodes(columns, master_barcode_file,
    #                             barcode_mapper_output)
    # elif master_barcode_file is not None:
    #     columns = unmap_barcodes(columns, master_barcode_file)

    df = pd.DataFrame(mat.todense(), columns=columns, index=rows)
    return df


def read_sv_as_dataframe(fn, sep, master_barcode_file, barcode_mapper_output):
    """
    Assumes tabbed file includes rownames (first column) and colnames (first row).

    :param fn: string
        filename
    :return:
    """
    df = pd.read_table(fn, index_col=0, sep=sep)

    # if master_barcode_file is not None and barcode_mapper_output is not None:
    #     df.columns = map_barcodes(df.columns, master_barcode_file,
    #                             barcode_mapper_output)
    # elif master_barcode_file is not None:
    #     df.columns = unmap_barcodes(df.columns, master_barcode_file)

    return df, df.columns, df.index


def write_dataframe_as_h5(df, output_file, genome, gene_ids, gene_names,
                          barcodes):
    h5_mtx = scipy.sparse.csc_matrix(df.values)


    flt = tables.Filters(complevel=1)
    with tables.open_file(output_file, 'w', filters=flt) as f:
        f.set_node_attr(f.root, "chemistry_description", "Single Cell 3\' V2")
        f.set_node_attr(f.root, "filetype", "matrix")
        # f.set_node_attr(f.root, "library_ids", ['id1'])
        # f.set_node_attr(f.root, "original_gem_groups", [1])
        f.root._f_setattr('original_gem_groups', np.array([1]))
        f.root._f_setattr('library_ids', np.array(['id1']))
        try:
            group = f.create_group(f.root, genome)
            # for attribute in ('indices', 'indptr'):
            #     arr = np.ndarray(getattr(h5_mtx, attribute))
            #     f.create_carray(group, attribute, obj=arr)
            f.create_carray(group, 'data', obj=np.asarray(h5_mtx.data, dtype=np.dtype('int32')))
            f.create_carray(group, 'genes', obj=np.asarray(gene_ids, dtype=np.dtype('S18')))
            f.create_carray(group, 'gene_names', obj=np.asarray(gene_names, dtype=np.dtype('S18')))
            f.create_carray(group, 'barcodes', obj=np.asarray(barcodes, dtype=np.dtype('S18')))
            f.create_carray(group, 'indices', obj=np.asarray(h5_mtx.indices, dtype=np.dtype('uint32')))
            f.create_carray(group, 'indptr', obj=np.asarray(h5_mtx.indptr, dtype=np.dtype('uint32')))
            f.create_carray(group, 'shape', obj=np.array(getattr(h5_mtx, 'shape'), dtype=np.dtype('int32')))
        except Exception as e:
            print('cannot write h5', e)


def write_dataframe_as_matrix(df, output_file):
    """
    Writes the pandas dataframe as a matrix.

    :param df:
    :param output_file:
    :return:
    """
    mtx = scipy.sparse.csr_matrix(df.values)
    scipy.io.mmwrite(output_file, mtx)


def write_dataframe_as_sv(df, output_file, sep='\t'):
    """
    Writes a {sep}-separated-values (sv) dataframe to file

    :param df:
    :param output_file:
    :param sep:
    :return:
    """
    df.to_csv(output_file, sep=sep)


def convert(input_file, input_type, output_file, output_type,
            columns_file, rows_file, genome,
            master_barcode_file, barcode_mapper_output):
    """
    Runs the conversion between input_type to output_type

    :param input_file:
    :param input_type:
    :param output_file:
    :param output_type:
    :param columns_file:
    :param rows_file:
    :return:
    """
    if input_type == 'csv':
        if output_type == 'mtx':
            print('warning: gene name+id may not be present in h5 file')
            sv2mtx(input_file, output_file, ',', master_barcode_file, barcode_mapper_output)
        elif output_type == 'h5':
            print('warning: gene name+id may not be present in h5 file')
            sv2h5(input_file, output_file, genome, ',', master_barcode_file, barcode_mapper_output)
    elif input_type == 'tsv':
        if output_type == 'mtx':
            print('warning: gene name+id may not be present in h5 file')
            sv2mtx(input_file, output_file, '\t', master_barcode_file, barcode_mapper_output)
        elif output_type == 'h5':
            print('warning: gene name+id may not be present in h5 file')
            sv2h5(input_file, output_file, genome, '\t', master_barcode_file, barcode_mapper_output)
    elif input_type == 'mtx':
        if output_type == 'tsv':
            mtx2sv(input_file, rows_file, columns_file, output_file, '\t', master_barcode_file, barcode_mapper_output)
        elif output_type == 'csv':
            mtx2sv(input_file, rows_file, columns_file, output_file, ',', master_barcode_file, barcode_mapper_output)
    elif input_type == 'h5':
        if output_type == 'tsv':
            h52sv(input_file, output_file, genome, '\t', master_barcode_file, barcode_mapper_output)
        elif output_type == 'csv':
            h52sv(input_file, output_file, genome, ',', master_barcode_file, barcode_mapper_output)


def sv2mtx(input_file, output_file, sep,
           master_barcode_file=None, barcode_mapper_output=None):
    """
    Converts a tab/comma/sep separated file into an mtx file

    :param input_file:
    :param output_file:
    :return:
    """
    df, barcodes, genes = read_sv_as_dataframe(input_file, sep, master_barcode_file, barcode_mapper_output)
    write_dataframe_as_matrix(df, output_file)

    o = open(output_file + '.columns', 'w')
    for column in barcodes:
        o.write(column + '\n')
    o.close()

    o = open(output_file + '.rows', 'w')
    for row in genes:
        o.write(row + '\n')
    o.close()


def mtx2sv(input_file, rows_file, columns_file, output_file, sep,
           master_barcode_file=None, barcode_mapper_output=None):
    """
    Converts mtx file + rows + columns into a tab/comma-separated file.

    :param input_file:
    :param rows_file:
    :param columns_file:
    :param output_file:
    :return:
    """
    df = read_mtx_as_dataframe(
        input_file, columns_file, rows_file, master_barcode_file, barcode_mapper_output
    )
    write_dataframe_as_sv(df, output_file, sep=sep)


def mtx2h5(input_file, output_file, columns_file, rows_file, genome,
           master_barcode_file=None, barcode_mapper_output=None):
    """
    Converts mtx/mtx associated files into h5 format.

    :param input_file:
    :param output_file:
    :param columns_file:
    :param rows_file:
    :param genome:
    :return:
    """
    df = read_mtx_as_dataframe(
        input_file, columns_file, rows_file, master_barcode_file, barcode_mapper_output
    )
    columns = [
        row[0] for row in csv.reader(open(columns_file), delimiter="\t")
        ]
    rows = [
        row[0] for row in csv.reader(open(rows_file), delimiter="\t")
        ]
    write_dataframe_as_h5(df, output_file, genome, rows, rows, columns)


def h52mtx(input_file, output_file, genome,
           master_barcode_file=None, barcode_mapper_output=None):
    """
    Converts from h5 file to mtx

    :param input_file:
    :param output_file:
    :param genome:
    :return:
    """
    gene_ids, gene_names, barcodes, matrix = get_matrix_from_h5(
        input_file, genome, master_barcode_file, barcode_mapper_output
    )

    if master_barcode_file is not None and barcode_mapper_output is not None:
        barcodes = map_barcodes(barcodes, master_barcode_file, barcode_mapper_output)
    elif master_barcode_file is not None:
        barcodes = unmap_barcodes(barcodes, master_barcode_file)

    o = open(output_file + '.rows', 'w')
    for i in range(0, len(gene_ids)):
        o.write('{}\t{}\n'.format(gene_ids[i], gene_names[i]))
    o.close()

    o = open(output_file + '.columns', 'w')
    for column in barcodes:
        o.write(column + '\n')
    o.close()

    mtx = scipy.sparse.csr_matrix(matrix)
    scipy.io.mmwrite(output_file, mtx)


def sv2h5(input_file, output_file, genome, sep,
          master_barcode_file=None, barcode_mapper_output=None):
    """
    Convert tab or csv separated files to h5.

    :param input_file:
    :param output_file:
    :param genome:
    :param sep:
    :return:
    """
    df, barcodes, gene_ids = read_sv_as_dataframe(input_file, sep, master_barcode_file, barcode_mapper_output)

    if master_barcode_file is not None and barcode_mapper_output is not None:
        barcodes = map_barcodes(barcodes, master_barcode_file, barcode_mapper_output)
    elif master_barcode_file is not None:
        barcodes = unmap_barcodes(barcodes, master_barcode_file)

    write_dataframe_as_h5(df, output_file, genome, gene_ids, gene_ids,
                          barcodes)


def h52sv(input_file, output_file, genome, sep,
          master_barcode_file=None, barcode_mapper_output=None):
    """
    Converts an h5 input file into a tab separated file.

    :param input_file:
    :param output_file:
    :param genome:
    :return:
    """
    gene_ids, gene_names, barcodes, mat = get_matrix_from_h5(input_file, genome, master_barcode_file, barcode_mapper_output)


    if master_barcode_file is not None and barcode_mapper_output is not None:
        barcodes = map_barcodes(barcodes, master_barcode_file, barcode_mapper_output)
    elif master_barcode_file is not None:
        barcodes = unmap_barcodes(barcodes, master_barcode_file)

    df = pd.DataFrame(mat.todense(), columns=barcodes, index=gene_ids)
    write_dataframe_as_sv(df, output_file, sep)


def map_barcodes(barcodes, master_barcode_file, barcode_mapper_output):
    """

    :param barcodes: list

    :param master_barcode_file: string
        file that contains the master barcode
    :param barcode_mapper_output: string
        file to write barcode mappings to
    :return:
    """
    # sorted_barcodes = sorted(barcodes)
    new_mapped_barcodes = []
    o = open(barcode_mapper_output, 'w')
    with open(master_barcode_file, 'r') as f:
        for current_barcode in barcodes: # sorted_barcodes:
            new_barcode = f.readline().rstrip()
            o.write('{}\t{}\n'.format(current_barcode, new_barcode))
            new_mapped_barcodes.append(new_barcode + '-1')
    o.close()
    return new_mapped_barcodes


def unmap_barcodes(barcodes, master_barcode_file):
    """
    Takes a barcodes master file and converts current set of barcodes back
    to the original

    :param barcodes:
    :param master_barcode_file:
    :return:
    """
    barcodes_dict = defaultdict()
    unmapped_barcodes_list = []

    with open(master_barcode_file, 'r') as f:
        for line in f:
            original_barcode, current_barcode = line.rstrip().split('\t')
            barcodes_dict[current_barcode] = original_barcode

    for barcode in barcodes:
        # print("old barcode: {}, new barcode: {}".format(barcode, barcodes_dict[barcode]))
        unmapped_barcodes_list.append(barcodes_dict[barcode])

    return unmapped_barcodes_list

def main():
    """
    Main program.

    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        required=True,
    )
    parser.add_argument(
        "--input_type",
        required=True,
        default='tsv'
    )
    parser.add_argument(
        "--output",
        required=True,
    )
    parser.add_argument(
        "--output_type",
        required=True,
        default='mtx'
    )
    parser.add_argument(
        "--columns",
        required=False,
        default=None,
        help="barcodes.tsv"
    )
    parser.add_argument(
        "--rows",
        required=False,
        default=None,
        help="genes.tsv"
    )
    parser.add_argument(
        "--genome",
        required=False,
        default=None,
        help="genome (hg19, mm10, etc.)"
    )
    parser.add_argument(
        "--master_barcode_file",
        required=False,
        default=None,
        help="master mapper"
    )
    parser.add_argument(
        "--barcode_mapper_output",
        required=False,
        default=None,
        help="master mapper output file"
    )

    args = parser.parse_args()

    input_file = args.input
    input_type = args.input_type
    output_file = args.output
    output_type = args.output_type
    columns_file = args.columns
    rows_file = args.rows
    genome = args.genome
    master_barcode_file = args.master_barcode_file
    barcode_mapper_output = args.barcode_mapper_output

    convert(
        input_file, input_type, output_file, output_type,
        columns_file, rows_file, genome,
        master_barcode_file, barcode_mapper_output
    )

if __name__ == "__main__":
    main()
