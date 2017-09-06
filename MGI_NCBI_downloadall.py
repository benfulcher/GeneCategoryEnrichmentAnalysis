#-------------------------------------------------------------------------------
# Idea is to download all MGI IDs listed for mouse, along with their
# NCBI (Entrez) Gene IDs (.csv output)
# (For info on Entrez IDs):
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3013746/
# cf. http://www.mousemine.org/mousemine/query.do#showing
#-------------------------------------------------------------------------------

from intermine.webservice import Service
service = Service("http://www.mousemine.org/mousemine/service")

# query description - Returns all genes for the specified organism.

# Get a new query on the class (table) you will be querying:
query = service.new_query("Gene")

# The view specifies the output columns
query.add_view("primaryIdentifier", "symbol", "ncbiGeneNumber")

# Uncomment and edit the line below (the default) to select a custom sort order:
query.add_sort_order("Gene.primaryIdentifier", "ASC")

# You can edit the constraint values below
query.add_constraint("organism.shortName", "=", "M. musculus", code = "A")
query.add_constraint("primaryIdentifier", "CONTAINS", "MGI:", code = "B")

# Uncomment and edit the code below to specify your own custom logic:
query.set_logic("A and B")

# Fill the MGI ID dictionary:
MGIDict = []
for row in query.rows():
    MGIDict.append({"MGIID":row["primaryIdentifier"], "symbol":row["symbol"], "NCBIGeneNumber":row["ncbiGeneNumber"]})
    # print row["primaryIdentifier"], row["symbol"], row["ncbiGeneNumber"]

# Convert to dataframe:
df = pd.DataFrame(MGIDict)

# Save out to csv
allDataFilename = "ALL_MGI_ID_NCBI.csv"
df.to_csv(allDataFilename)
