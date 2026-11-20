from harvest_ncbi_core_metadata import *
# Identified from the references of the ATCC20541 comparative-genome paper.
pmids=['31806745','32299870','34485249','30394102']
eu('efetch',{'db':'pubmed','id':','.join(pmids),'retmode':'xml'},'reference_followup_pubmed.xml','xml')
eu('efetch',{'db':'pubmed','id':'32108289,41497419','retmode':'xml'},'verified_followup_pubmed.xml','xml')
eu('efetch',{'db':'pmc','id':'12766970','retmode':'xml'},'TT0037_publication.xml','xml')
eu('efetch',{'db':'pmc','id':'6895305,7163008,8416056','retmode':'xml'},'reference_followup_pmc.xml','xml')
print('Additional citation metadata downloaded',flush=True)
