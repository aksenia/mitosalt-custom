"""
SaltShaker CLI

Command-line interface with subcommands for modular analysis.
"""

import argparse
import sys
from .cli import call, classify, plot


def main():
    parser = argparse.ArgumentParser(
        prog='saltshaker',
        description='SaltShaker: Pattern classification and visualization for MitoSAlt'
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Add subcommand parsers
    call.add_parser(subparsers)
    classify.add_parser(subparsers)
    plot.add_parser(subparsers)
    
    # Parse arguments
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        sys.exit(1)
    
    # Execute the appropriate subcommand
    if args.command == 'call':
        call.run(args)
    elif args.command == 'classify':
        classify.run(args)
    elif args.command == 'plot':
        plot.run(args)


if __name__ == "__main__":
    main()