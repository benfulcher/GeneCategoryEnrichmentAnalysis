from intermine.webservice import Service
import pandas as pd
import csv # For saving string data to csv

# Idea is to convert a list of MGI IDs (.csv input) -> NCBI (Entrez) Gene IDs (.csv output)
# (For info on Entrez IDs):
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3013746/

# Initiate a mousemine service:
service = Service("http://www.mousemine.org/mousemine/service")

# Read in MGI IDs:
MGIIDFile = "MGI_IDs.csv"
with open(MGIIDFile,'rb') as f:
    reader = csv.reader(f)
    MGIIDList = list(reader)
MGIIDList = [MGIIDList[x][0] for x in range(len(MGIIDList))]

# Convert to NCBI Gene IDs
MGIDict = []
for x in range(len(MGIIDList)):
    query = service.new_query("Gene")
    query.add_view("primaryIdentifier", "ncbiGeneNumber", "symbol")
    query.add_constraint("primaryIdentifier", "=", MGIIDList[x], code = "A")
    for row in query.rows():
        print x, '/', len(MGIIDList), ':', row["primaryIdentifier"], row["symbol"], row["ncbiGeneNumber"]
        MGIDict.append({'MGIID':MGIIDList[x], 'symbol':row['symbol'], 'NCBIGeneNumber':row["ncbiGeneNumber"]})
    if x>0 and ((x % 500)==0):
        # To dataframe:
        df = pd.DataFrame(MGIDict)

        # Save out:
        allDataFilename = "MGI_ID_NCBI_%u.csv" % x
        df.to_csv(allDataFilename)
        print 'Saved yo - %u' % x

# To dataframe:
df = pd.DataFrame(MGIDict)

# Save out:
allDataFilename = "MGI_ID_NCBI.csv"
df.to_csv(allDataFilename)
