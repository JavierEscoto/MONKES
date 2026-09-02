
#!/usr/bin/env python3

import sys
from read_configuration import read_configuration
from write_header import write_header
from create_monkes_runs import create_monkes_runs
from read_monkes_database import read_monkes_database

def main():
    if len(sys.argv) < 2:
        print("Usage: create_dk_database.py <filename>")
        sys.exit(1)

    filename = sys.argv[1]

    # Step 1: Read equilibrium information from the provided filename
    params = read_configuration(filename)

    # Step 2: Write monkes.dk header
    write_header(params)

    # Step 3: Create MONKES input directories/files and submit jobs
    create_monkes_runs(params,filename)

    # Step 4: Read MONKES output once all runs complete
    read_monkes_database(params)

if __name__ == "__main__":
    main()

