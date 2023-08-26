import subprocess
'''
This function is used to automitize the process of data extraction from the VCF files of the Ensembl Variation database, using the
bcftools program.
'''

# Define the range of numbers you want to iterate over
for i in range(3, 22):
    # Construct the output filename with the current number
    output_filename = f"/mnt/sda1/Databases/Ensembl/Variation/110/chromosomes_data/chr{i}_data.tsv"

    # Construct the input VCF filename with the current number
    input_filename = f"/mnt/sda1/Databases/Ensembl/Variation/110/VCF/homo_sapiens-chr{i}.vcf.gz"

    # Construct the bcftools query command
    cmd = [
        "bcftools",
        "query",
        "-f", '%CHROM\t%POS\t%REF\t%ALT\t%TSA\t%ID',
        "-H",
        "-o", output_filename,
        input_filename
    ]

    # Run the command
    subprocess.run(cmd)
'''
# Construct the output filename with the current number
output_filename = f"/mnt/sda1/Databases/Ensembl/Variation/110/chromosomes_data/chr22_data.tsv"

# Construct the input VCF filename with the current number
input_filename = f"/mnt/sda1/Databases/Ensembl/Variation/110/VCF/homo_sapiens-chr22.vcf.gz"

# Construct the bcftools query command
cmd = [
    "bcftools",
    "query",
    "-f", "'%CHROM\t%POS\t%REF\t%ALT\t%TSA\t%ID'",
    "-H",
    "-o", output_filename,
    input_filename
]

# Run the command
subprocess.run(cmd)
'''

print("Processing complete.")