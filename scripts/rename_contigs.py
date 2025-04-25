import argparse

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', '--input', help='input FASTA', required=True)
  parser.add_argument('-b', '--bedfile', help='bedfile with contig', required=False)
  parser.add_argument('-o', '--output', help='output FASTA', required=True)
  parser.add_argument('-c', '--contigs', help='contig names', required=False, nargs='+')
  args = parser.parse_args()

  new_names = []

  if args.bedfile:
    bedfile = open(args.bedfile)
    bedlines = bedfile.readlines()
    # Make sure bedfiles end with extra line
    for i in bedlines:
      new_names.append(i.split('	')[3])
    bedfile.close()
  elif args.contigs:
    new_names.append(args.contigs)
    new_names=new_names[0]
    new_names=[i + "\n" for i in new_names]
  else:
    print("Either a bedfile or array of contig names must be specified")
    exit()

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