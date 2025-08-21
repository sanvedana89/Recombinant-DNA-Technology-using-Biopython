from Bio.Restriction import *
from Bio.Restriction.Restriction import RestrictionType

# List all enzyme names
enzyme_list = [enzyme for enzyme in AllEnzymes]

# Display details for each enzyme
for enzyme in enzyme_list[:1000]:  # Display first 10 enzymes (to avoid too much output)
    print(f"Enzyme: {enzyme.__name__}")
    print(f"   Recognition site: {enzyme.site}")
    print(f"   Overhang: {enzyme.ovhg}")
    print(f"   Isoschizomers: {enzyme.isoschizomers()}")
    print("-" *10 )
