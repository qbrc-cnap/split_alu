import argparse
import regex

# save the write operation until we have this many fastq records:
BUFFER_SIZE = 10000

def get_rc(s):
    mapping = dict(zip(list('ACGT'), list('TGCA')))
    rc = ''.join([mapping[x] for x in s[::-1]])
    return rc


def read_cl():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Input FASTQ file', required=True)
    parser.add_argument('-o', '--output', help='Output FASTQ file', required=True)
    parser.add_argument('-b', '--barcode', help='Barcode/tag sequence', required=True)
    args = parser.parse_args()
    return vars(args)


def search_seq(s, barcode, rc):
    i = 0
    locations = [0,]
    for x in regex.finditer('(?e)(%s)|(%s){e<=1}' % (barcode,rc), s):
        match_start = x.start()
        match_end = x.end()
        locations.append(x.start())
        locations.append(x.end())
    locations.append(len(s))
    return locations


def get_records(fq):
    fastq_record = None
    for i, line in enumerate(open(fq)):
        if i % 4 == 0:
            yield fastq_record
            fastq_record = []
            readname = line.strip().split(' ')[0]
            fastq_record.append(readname)
        else:
            fastq_record.append(line.strip())
    yield fastq_record


if __name__ == '__main__':
    args = read_cl()
    barcode = args['barcode']
    rc = get_rc(barcode)
    with open(args['output'], 'w') as fout:
        buffer_count= 0
        lines = []
        for fastq_record in get_records(args['input']):
            if fastq_record:
                breakpoint_list = search_seq(fastq_record[1], barcode, rc)
                for i in range(0, len(breakpoint_list), 2):
                    a = breakpoint_list[i]
                    b = breakpoint_list[i+1]
                    seq = fastq_record[1][a:b]
                    qual = fastq_record[3][a:b]
                    readname = '%s?%d/%d:%d-%d' % (fastq_record[0], i//2, len(breakpoint_list)//2-1 , a, b)
                    if b > a:
                        lines.extend([readname, seq, '+', qual])
                    buffer_count += 1
            if (buffer_count % BUFFER_SIZE) == 0:
                if len(lines) > 0:
                    fout.write('\n'.join(lines) + '\n')
                    lines = []
                    buffer_count = 0
        fout.write('\n'.join(lines) + '\n')
