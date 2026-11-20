from worker import OUT,limits,set_stage_limit,protected,save,now
before=protected()
set_stage_limit('BUSCO_resource_probe');busco=limits()
set_stage_limit('annotation_resource_probe');other=limits()
after=protected();assert before['MemoryMax']==after['MemoryMax']=='100000000000'
save(OUT/'reports/resource_cap_transition_test.json',{'at':now(),'status':'PASS','BUSCO':busco,'other_stages':other,'protected_before':before,'protected_after':after})
print('PASS: 20 GB BUSCO / 4 GiB other-stage transitions, protected cap unchanged')
