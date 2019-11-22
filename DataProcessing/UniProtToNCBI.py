import urllib,urllib2
import csv # For saving string data to csv
import pandas as pd
url = 'http://www.uniprot.org/uploadlists/'

#-------------------------------------------------------------------------------
# Read in UniProt IDs from file:
UniProtIDFile = "allUniprotIDs.csv"
with open(UniProtIDFile,'rb') as f:
    reader = csv.reader(f)
    UniProtIDList = list(reader)
UniProtIDList = [UniProtIDList[x][0] for x in range(len(UniProtIDList))]

#-------------------------------------------------------------------------------
UniProtDict = []
for x in range(len(UniProtIDList)):
    params = {
    'from':'ID',
    'to':'P_ENTREZGENEID',
    'format':'list',
    'query':UniProtIDList[x]
    }

    data = urllib.urlencode(params)
    request = urllib2.Request(url, data)
    contact = "" # Please set your email address here to help us debug in case of problems.
    request.add_header('User-Agent', 'Python %s' % contact)
    response = urllib2.urlopen(request)
    page = response.read(200000)
    page = page.replace("\n", " ").strip()
    print x, '/', len(UniProtIDList), ':', UniProtIDList[x], '->', page
    UniProtDict.append({'UniProtID':UniProtIDList[x], 'EntrezID':page})

df = pd.DataFrame(UniProtDict)

# Save out:
allDataFilename = "UniProt_Entrez_Map.csv"
df.to_csv(allDataFilename)
