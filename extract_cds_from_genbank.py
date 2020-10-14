#!/usr/bin/env python3
from pathlib import Path

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import SeqRecord

gb_file = Path("AE016877.1.gb")


def get_genbank(query: str, retmax=10):
    # find these values with einfo
    Entrez.email = "Yi.Zhou@uga.edu"

    # Entrez esearch: search query in specified database
    r = Entrez.read(Entrez.esearch(db="assembly", term=query, retmax=retmax))
    taxid = r["IdList"]

    # Entrez elink: link assembly accession number to nucleotide accession
    r = Entrez.read(Entrez.elink(from_uid=taxid, dbfrom="assembly", db="nucleotide"))
    nuccore_insdc = [
        x for x in r[0]["LinkSetDb"] if x["LinkName"] == "assembly_nuccore_insdc"
    ]
    dl_links = [x["Id"] for x in nuccore_insdc[0]["Link"]]

    # Entrez efetch: download GenBank sequences
    records = []
    for dl_link in dl_links:
        handle = Entrez.efetch(
            id=dl_link, db="nucleotide", rettype="gbwithparts", retmode="text"
        )
        records.append(SeqIO.read(handle, "genbank"))

    return records


if not gb_file.is_file():
    gb = get_genbank(query="AE016877.1", retmax=1)
else:
    gb = SeqIO.parse(gb_file.absolute(), "genbank")

genome = [x for x in gb if x.name == "AE016877"][0]  # drop the plasmid
features = [
    x for x in genome.features if "locus_tag" in x.qualifiers and x.type == "CDS"
]  # throw away gene, operon, rRNA, ...

# Locus tags of the four genes in the operon we're interested in
operon_tags = ["BC_3514", "BC_3515", "BC_3516", "BC_3517"]
operon = [x for x in features if x.qualifiers["locus_tag"][0] in operon_tags]

# We want to study 10kb flanking the operon
start_pos: int = min([min(x.location) for x in operon])
end_pos: int = max([max(x.location) for x in operon])

operon_range = FeatureLocation(start=start_pos - 10000, end=end_pos + 10000)

# Get all the proteins in the flanking region (including the operon itself)
res = []
for f in features:
    if f.location.start in operon_range or f.location.end in operon_range:
        res.append(f)

# Prepare output as fasta file
out_file = []
for record in res:
    protein_seq = record.qualifiers["translation"][0]
    record_id = record.qualifiers["locus_tag"][0]
    record_name = record.qualifiers["protein_id"][0]
    record_description = record.qualifiers["product"][0]
    out_file.append(
        SeqRecord(
            seq=Seq(protein_seq),
            id=record_id,
            name=record_name,
            description=record_description,
        )
    )

SeqIO.write(out_file, "AE016877.1.faa", "fasta")
