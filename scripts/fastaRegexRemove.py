import argparse
import re

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', '--input', help='input FASTA', required=True)
  parser.add_argument('-r', '--regex', help='regex string in names to remove', required=True)
  parser.add_argument('-o', '--output', help='output FASTA', required=True)
  parser.add_argument('-n', '--new_string', help='string to replace regex', required=False, default='', type=str)
  args = parser.parse_args()

  in_fasta = open(args.input)
  out_fasta = open(args.output, 'w')

  args.new_string if args.new_string != 'default_val' else ''

  in_fasta_lines = in_fasta.readlines()
  for i in in_fasta_lines:
    if i.startswith('>'):
      contig_name = re.sub(args.regex, args.new_string, i)
      out_fasta.write(contig_name)
    else:
      out_fasta.write(i)
  
  out_fasta.close()
  in_fasta.close()

if __name__ == '__main__':
  main()
  