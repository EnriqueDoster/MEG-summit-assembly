#!/usr/bin/env python3

import gdown
import sys
import argparse

def parse_cmdline_params(cmdline_params):
    info = ""
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-u', '--url', nargs='+', required=True,
                        help='URL to file in google drive')
    parser.add_argument('-o', '--output_file', required=True,
                        help='Output directory for writing the AMR_analytic_matrix.csv file')
    return parser.parse_args(cmdline_params)

if __name__ == '__main__':
    opts = parse_cmdline_params(sys.argv[1:])
    url = opts.url
    print(url)
    output = opts.output_file
    print(output)
    gdown.download(url,output,quiet=False)
    #gdown.download(url,output,quiet=False)
