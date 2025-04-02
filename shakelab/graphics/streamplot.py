import argparse
import sys
import os
import glob

VERSION = '202400822'


def streamplot(input_paths, file_format):
    """
    Process the input files or directories and generate plots.

    Parameters:
    input_paths (list): A list of paths to the input files or directories.
    file_format (str): The format of the input files.
    """
    all_files = []
    for path in input_paths:
        if os.path.isdir(path):
            # If the path is a directory, gather all files in the directory
            search_pattern = os.path.join(path, "*")
            matched_files = glob.glob(search_pattern)
            if not matched_files:
                print(f"Warning: No files found in directory {path}")
            else:
                all_files.extend(matched_files)
        else:
            # If the path is a file or pattern, gather matching files
            matched_files = glob.glob(path)
            if not matched_files:
                print(f"Warning: No files matched the pattern {path}")
            else:
                all_files.extend(matched_files)
    
    if not all_files:
        print("Error: No valid input files found.")
        sys.exit(1)

    for input_file in all_files:
        if not os.path.isfile(input_file):
            print(f"Error: File not found: {input_file}")
            continue
        print(f"Processing file: {input_file} with format: {file_format}")
        # Add your file processing and plotting code here
    pass


if __name__ == '__main__':
    info_string = 'Plot utility for seismic data'

    parser = argparse.ArgumentParser(
        description=info_string,
        usage='python %(prog)s [options] input_files_or_directories...',
        epilog=('Example: python script.py -f ascii data_directory '
                'file1.ascii file2.ascii')
    )

    # Adding an option to display the version
    parser.add_argument('-v', '--version', action='version',
                        version=f'Version {VERSION}')

    # Adding the argument for the input files or directories
    parser.add_argument('input_paths', nargs='+', metavar='input_paths',
                        type=str, help=('The input files or directories '
                                        '(supports wildcards)'))

    # Adding the option to specify the file format
    parser.add_argument('-f', '--file_format', type=str, default='',
                        help='Specify file format (optional)')

    # Parsing the arguments
    args = parser.parse_args()

    # Call the streamplot function with the arguments
    streamplot(args.input_paths, args.file_format)