from intermine.webservice import Service
import pandas as pd
import csv # For saving string data to csv

def SaveListCSV(stringList,fileName):
    # Outputs a csv from a given list of strings
    resultFile = open(fileName,'wb')
    wr = csv.writer(resultFile, dialect='excel')
    wr.writerow(stringList)


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

# To dataframe:
df = pd.DataFrame(MGIDict)

# Save out:
allDataFilename = "MGI_ID_NCBI.csv"
df.to_csv(allDataFilename)
