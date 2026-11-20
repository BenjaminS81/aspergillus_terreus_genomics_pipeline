import annotate_terreus_trees_v2 as p
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
fig,ax=plt.subplots(figsize=(13,10));ax.axis('off');ax.set_title('Shared species colours and annotation policy',loc='left',fontsize=19,weight='bold',color='#213b54')
ax.legend(handles=[Patch(color=c,label=s) for s,c in p.colors.items()],loc='upper left',ncol=2,frameon=False,fontsize=10)
fig.text(.07,.075,'Species ring: curated reference labels, documented ex-type genome cultures, or three-locus supported project placement.\nUnresolved project labels remain explicit; no unsupported sensu-stricto or formal subspecies assignment.\n[T] External ex-type reference. [TG] Documented genome ex-type culture; see the provenance audit.\nReference type-clade ring: SH-aLRT >=80 and UFBoot >=95. Genome-type evidence has its own outer ring.\nType-focused figures prune the full trees; genome-only reference rings are projected and explicitly labelled.',fontsize=10,color='#526778')
p.save(fig,'annotation_legend')
