from intermine.webservice import Service
service = Service("http://www.mousemine.org/mousemine/service")
query = service.new_query("Gene")
query.add_view("primaryIdentifier", "ncbiGeneNumber", "symbol")
query.add_constraint("primaryIdentifier", "=", "MGI:1918911", code = "A")

for row in query.rows():
    print row["primaryIdentifier"], row["symbol"], row["ncbiGeneNumber"]
