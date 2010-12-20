#!/usr/bin/env python
import common
import table
import hmm
import data

def main():
    tab = table.load_seq_data()
    h = hmm.init_hmm(tab, 5)
    for row in tab:
        print h.call_row(row)

if __name__ == '__main__':
    main()
