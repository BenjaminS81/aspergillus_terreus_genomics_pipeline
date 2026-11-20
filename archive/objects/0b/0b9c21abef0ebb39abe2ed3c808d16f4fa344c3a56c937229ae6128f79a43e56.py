"""Shared, resource-bounded helpers for the November 13 figure review."""
from pathlib import Path
import os, json, ctypes
BASE = Path(__file__).resolve().parent
if BASE.name=='scripts' and (BASE.parent/'data').is_dir():
    OUT=BASE.parent;OLD=OUT.parent/'2026-11-11'
else:
    OLD = BASE / 'comparative_visualization/2026-11-11'
    OUT = BASE / 'comparative_visualization/2026-11-13'

def constrain_process():
    """Limit this Windows worker only; never inspect or alter another job."""
    for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
        os.environ[key] = '1'
    if os.name != 'nt':
        return {'threads':1,'memory_cap':'not applied on this platform'}
    from ctypes import wintypes as w
    class BASIC(ctypes.Structure):
        _fields_=[('PerProcessUserTimeLimit',ctypes.c_int64),('PerJobUserTimeLimit',ctypes.c_int64),('LimitFlags',w.DWORD),('MinimumWorkingSetSize',ctypes.c_size_t),('MaximumWorkingSetSize',ctypes.c_size_t),('ActiveProcessLimit',w.DWORD),('Affinity',ctypes.c_size_t),('PriorityClass',w.DWORD),('SchedulingClass',w.DWORD)]
    class IO(ctypes.Structure):
        _fields_=[(k,ctypes.c_uint64) for k in ('ReadOperationCount','WriteOperationCount','OtherOperationCount','ReadTransferCount','WriteTransferCount','OtherTransferCount')]
    class EXT(ctypes.Structure):
        _fields_=[('BasicLimitInformation',BASIC),('IoInfo',IO),('ProcessMemoryLimit',ctypes.c_size_t),('JobMemoryLimit',ctypes.c_size_t),('PeakProcessMemoryUsed',ctypes.c_size_t),('PeakJobMemoryUsed',ctypes.c_size_t)]
    kernel=ctypes.WinDLL('kernel32',use_last_error=True)
    kernel.CreateJobObjectW.argtypes=[ctypes.c_void_p,w.LPCWSTR];kernel.CreateJobObjectW.restype=w.HANDLE
    kernel.GetCurrentProcess.restype=w.HANDLE
    kernel.SetInformationJobObject.argtypes=[w.HANDLE,ctypes.c_int,ctypes.c_void_p,w.DWORD]
    kernel.AssignProcessToJobObject.argtypes=[w.HANDLE,w.HANDLE]
    kernel.GetProcessAffinityMask.argtypes=[w.HANDLE,ctypes.POINTER(ctypes.c_size_t),ctypes.POINTER(ctypes.c_size_t)]
    kernel.SetProcessAffinityMask.argtypes=[w.HANDLE,ctypes.c_size_t]
    job=kernel.CreateJobObjectW(None,None); limits=EXT()
    limits.BasicLimitInformation.LimitFlags=0x200;limits.JobMemoryLimit=2*1024**3
    if not kernel.SetInformationJobObject(job,9,ctypes.byref(limits),ctypes.sizeof(limits)):
        raise ctypes.WinError(ctypes.get_last_error())
    process=kernel.GetCurrentProcess()
    if not kernel.AssignProcessToJobObject(job,process):raise ctypes.WinError(ctypes.get_last_error())
    own=ctypes.c_size_t();system=ctypes.c_size_t()
    if not kernel.GetProcessAffinityMask(process,ctypes.byref(own),ctypes.byref(system)):raise ctypes.WinError(ctypes.get_last_error())
    mask=own.value & -own.value
    if not kernel.SetProcessAffinityMask(process,mask):raise ctypes.WinError(ctypes.get_last_error())
    globals()['_job_handle']=job
    return {'job_memory_limit_bytes':2*1024**3,'cpu_affinity_mask':mask,'threads':1,'scope':'this worker only'}

def readjson(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def writejson(path,value):Path(path).write_text(json.dumps(value,indent=2),encoding='utf-8')
def fasta(path):
    name=None;chunks=[]
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith('>'):
                if name is not None:yield name,''.join(chunks).upper()
                name=line[1:].split()[0];chunks=[]
            else:chunks.append(line.strip())
    if name is not None:yield name,''.join(chunks).upper()

