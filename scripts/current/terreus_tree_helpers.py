import re
class Node:
    def __init__(self,name='',length=0,children=None):self.name=name;self.length=length;self.children=children or []
    def leaves(self):return [self] if not self.children else [x for c in self.children for x in c.leaves()]
def parse_newick(text):
    tokens=re.findall(r"[(),:;]|[^(),:;\s]+",text);pos=0
    def node():
        nonlocal pos
        n=Node()
        if tokens[pos]=='(':
            pos+=1;n.children.append(node())
            while tokens[pos]==',':pos+=1;n.children.append(node())
            assert tokens[pos]==')';pos+=1
        if tokens[pos] not in [',',')',':',';']:n.name=tokens[pos].strip("'");pos+=1
        if tokens[pos]==':':pos+=1;n.length=float(tokens[pos]);pos+=1
        return n
    result=node();assert tokens[pos]==';';return result
def support(n):
    try:
        v=[float(x) for x in n.name.split('/')]
        return v if len(v)==2 else [None,v[0]]
    except (ValueError,IndexError):return [None,None]
def strong(n):
    a,b=support(n);return a is not None and a>=80 and b>=95
def prune(n,keep):
    if not n.children:return Node(n.name,n.length) if n.name in keep else None
    children=[x for c in n.children if (x:=prune(c,keep)) is not None]
    if not children:return None
    if len(children)==1:
        c=children[0];c.length+=n.length
        if c.children:c.name=''
        return c
    return Node(n.name,n.length,children)
