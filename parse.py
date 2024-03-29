#!/usr/bin/env python3
import matplotlib.pyplot as plt
import argparse, tables
import pandas as pd
import numpy as np

def main():
    parser = argparse.ArgumentParser(description='Parses h5 files generated by monteCarlo and makes plots')
    parser.add_argument('-f', '--file', type=str, required=True, help='Input file')
    args = parser.parse_args()

    print(f'Reading {args.file}...', end='')
    infile = tables.open_file(args.file)
    params = infile.root.table._v_attrs # get run attributes
    df = pd.DataFrame.from_records( infile.root.table.read() ) # use either read() or read_where()
    infile.close()
    print('done. File closed')

    for attr in params._f_list("user"):
        print(f"{attr}: {getattr(params, attr)}")

    df['status'] = df['status'].str.decode("utf-8") # Status column in byte string format

    print('Total neutrons simulated: ', len(df.index) )
    print('Neutrons lost at source: ', len( df.query('status=="source"').index ) )
    print('Neutrons lost in pipe: ', len( df.query('status=="pipe"').index ) )
    print('Neutrons lost on window: ', len( df.query('status=="window"').index ) )

    print("\n### Neutrons in the cell ###")
    df.query('status=="cell"', inplace=True)
    print("Total neutrons in cell: ", len(df.index))
    print("Av window hits = ", df['windowHits'].mean() )
    print("Av cell rejections = ", df['cellRejections'].mean() )
    print("Av total bounces= ", df['totalSteps'].mean() )

    plt.figure()
    plt.title('windowHits (end = cell)')
    plt.xlabel('Number of window hits')
    plt.ylabel('Number of neutrons')
    df['windowHits'].hist(bins= int(df['windowHits'].max() + 2 ), range=[0, int(df['windowHits'].max()) + 2], align='left' )
    plt.grid(False)
    plt.show()

    return


if ( __name__ == '__main__' ):
    main()
