import prodigal
import hmmer
from argparse import ArgumentParser
import os.path


def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg

def is_valid_folder(parser, arg):
    if not os.path.isdir(arg):
        parser.error("The folder %s does not exist!" % arg)
    else:
        return arg



def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("-E", "--evalue", help="evalue threshold",
                    type=double, default=1e-8)

    parser.add_argument("-i", "--input", dest="input_file", required=True,
                    help="input file",
                    type=lambda x: is_valid_file(parser, x))

    parser.add_argument("-p", "--train_dir", dest="train_dir",
                    help=("location of a training directory (a previous run" +
                        " of squeaky)"),
                    type=lambda x: is_valid_folder(parser, x))

    parser.add_argument("-o", "--out_dir", dest="output_dir", required=True,
                    help="location of an output directory",
                    type=lambda x: is_valid_folder(parser, x))

    args = parser.parse_args()



    return


if __name__ == '__main__':
    main()
