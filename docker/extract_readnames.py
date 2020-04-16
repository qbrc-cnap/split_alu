import argparse
import gzip


def read_cl():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r1', help='R1 fastq.gz file', required=True)
    parser.add_argument('-r2', help='R2 fastq.gz file', required=True)
    parser.add_argument('-o', '--output', help='output file', required=True)
    args = parser.parse_args()
    return vars(args)


def extract_names(r1, r2, output):
    s1 = set()
    s2 = set()
    for i, line in enumerate(gzip.open(r1)):
        if i % 4 == 0:
            l = line.decode('utf-8').strip().split(' ')[0][1:]
            s1.add(l)

    for i, line in enumerate(gzip.open(r2)):
        if i % 4 == 0:
            l = line.decode('utf-8').strip().split(' ')[0][1:]
            s2.add(l)

    u = s1.union(s2)
    fout = open(output, 'w')
    for x in u:
        fout.write('%s\n' % x)
    fout.close()


if __name__ == '__main__':
    args = read_cl()
    extract_names(args['r1'], args['r2'], args['output'])

