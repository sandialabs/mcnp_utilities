#!/usr/bin/env python3

from pypdf import PdfReader, PdfWriter
from argparse import ArgumentParser


def get_arguments():
  parser = ArgumentParser(description='Crop pdf produced by MCNP to just the figure')
  parser.add_argument(
    'pdf',
    type=str,
    help='Input PDF to crop',
    metavar='<PDF path>'
  )
  return parser.parse_args()

if __name__ == '__main__':
  args = get_arguments()
  reader, writer = PdfReader(args.pdf), PdfWriter()

  for i in range(len(reader.pages)):
    page = reader.pages[i]
    page.cropbox.lower_left = (81, 262)
    page.cropbox.upper_right = (531, 710)
    writer.add_page(page)

  [page.compress_content_streams(level=9) for page in writer.pages]

  with open(args.pdf, "wb") as fp:
    writer.write(fp)
