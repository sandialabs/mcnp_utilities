#! /usr/bin/env python3

from argparse import ArgumentParser


def convertKelvin(T):
  return 8.617E-11 * T

def convertCelsius(T):
  return convertKelvin(T + 273.15)

def convertFahrenheit(T):
  return 4.787E-11 * T

def convertRankine(T):
  return convertFahrenheit(T + 459.67)

def getConverter(Tunit):
  if Tunit == 'K':
    return convertKelvin
  elif Tunit == 'C':
    return convertCelsius
  elif Tunit == 'F':
    return convertFahrenheit
  elif Tunit == 'R':
    return convertRankine
  else:
    raise ValueError(f'Unrecognized temperature unit "{Tunit}"!')

def getArguments():
  parser = ArgumentParser(description="A tool to compute the TMP entry in MCNP for a given temperature")
  parser.add_argument(
    'temp',
    help='temperature value',
    type=float,
    metavar='<temp>'
  )
  parser.add_argument(
    'unit',
    help='Temperature unit (K=Kelvin, C=Celsius, F=Fahrenheit, R=Rankine)',
    type=str,
    choices=['K', 'C', 'F', 'R'],
    metavar='<unit>'
  )
  return parser,parser.parse_args()

if __name__ == '__main__':
  argp, args = getArguments()
  print(f'TMP={getConverter(args.unit)(args.temp):.5E}')
