from terreus_reference_sources import *
for uid in ['3161753','4260807','11736086','10142542','2681604']:
    r=eu('efetch',{'db':'pmc','id':uid,'retmode':'xml'},'PMC'+uid+'.xml')
    art=r.find('.//article');print('PAPER',uid,''.join(art.find('./front/article-meta/title-group/article-title').itertext()),'bytes',len(ET.tostring(r)),flush=True)
get('https://www.studiesinmycology.org/images/Volume/Sim107/Vol107Art1.pdf','Visagie_2024_review.pdf','bytes')
get('https://www.studiesinmycology.org/images/Volume/Sim107/Vol107Art1supp.zip','Visagie_2024_supplements.zip','bytes')
get('https://fuse-journal.org/images/Issues/Vol14Art6.pdf','Brazil_caves_2024.pdf','bytes')
eu('efetch',{'db':'pmc','id':'12786731','retmode':'xml'},'PMC12786731.xml')
print('PDF DOWNLOADED',flush=True)
