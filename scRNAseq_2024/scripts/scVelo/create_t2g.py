import sys

# Input GTF file
infile = sys.argv[1]
outfile = "t2g.txt"

# Process the input file
with open(infile, "r") as fin, open(outfile, "w") as fout:
    for line in fin:
        # Skip empty lines and comments
        if not line.strip() or line.startswith("#"):
            continue

        # Split the line into fields
        fields = line.strip().split("\t")
        
        # Ensure there are enough fields (e.g., 9 fields in a valid GTF line)
        if len(fields) < 9:
            continue
        
        # Check if the feature is a transcript
        if fields[2] != "transcript":
            continue

        # Extract attributes (the last column)
        attributes = fields[-1].strip()
        
        # Initialize default values
        gene_id = "NoGeneID"
        gene_name = ""

        # Parse attributes
        for attr in attributes.split(";"):
            attr = attr.strip()
            if attr.startswith("gene_id"):
                gene_id = attr.split()[1].strip('"')
            elif attr.startswith("gene_name"):
                gene_name = attr.split()[1].strip('"')

        # Write output: gene_id in the first column, gene_name or gene_id in the second column
        fout.write(f"{gene_id}\t{gene_name if gene_name else gene_id}\n")
