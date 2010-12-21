#!/usr/bin/env python
import common
import table
import hmm
import data

def main():
    tab = table.load_seq_data()
    h = hmm.train_hmm(tab, 5)

if __name__ == '__main__':
    main()
