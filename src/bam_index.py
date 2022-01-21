import sys
from lib import iter_bgzf_offsets


def main():
    with open(sys.argv[1], "rb") as bam, open(sys.argv[2], "wt") as out:
        for offset in iter_bgzf_offsets(bam):
            out.write(str(offset))
            out.write("\n")


if __name__ == "__main__":
    main()
