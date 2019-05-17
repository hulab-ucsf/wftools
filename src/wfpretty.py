import argparse
from bs4 import BeautifulSoup

g_version = "0.1"
g_exename = "wfpretty"


def getArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--filename", help="Input Filename", required=True)
    parser.add_argument("-o", "--output", help="Output Filename", required=True)
    return parser.parse_args()


def runApp(fn, outputFn):
    with open(fn, "r") as f:
        soup = BeautifulSoup(f, "xml")
        g = open(outputFn, "w")
        g.write(soup.prettify(encoding="utf-8").decode("utf-8"))
        g.close()

print("{0} v{1} - Copyright(c) HuLab@UCSF 2019".format(g_exename, g_version))
args = getArgs()
runApp(args.filename, args.output)
print("done.")
