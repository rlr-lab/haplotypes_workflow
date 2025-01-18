import argparse

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', '--input', help='input FASTA', required=True)
  parser.add_argument('-b', '--bedfile', help='bedfile with contig', required=True)
  parser.add_argument('-o', '--output', help='output FASTA', required=True)
  args = parser.parse_args()

  bedfile = open(args.bedfile)
  bedlines = bedfile.readlines()
  new_names = []
  for i in bedlines:
    new_names.append(i.split('	')[3])
  bedfile.close()

  in_fasta = open(args.input)
  out_fasta = open(args.output, 'w')

  name_count = 0
  in_fasta_lines = in_fasta.readlines()
  for i in in_fasta_lines:
    if i.startswith('>'):
      contig_name = '>' + new_names[name_count]
      out_fasta.write(contig_name)
      name_count += 1
    else:
      out_fasta.write(i)
  
  out_fasta.close()
  in_fasta.close()

if __name__ == '__main__':
  main()