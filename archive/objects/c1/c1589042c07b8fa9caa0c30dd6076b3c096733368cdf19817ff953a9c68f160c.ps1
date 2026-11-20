param([switch]$Resume)
# Default is verification only. Pass -Resume only when ready to resume analysis.
$resumeMode = if ($Resume) { '--resume' } else { '--check' }
wsl.exe -d Ubuntu -- /home/benjamin/.local/share/Achromobacter_ANI/envs/busco/bin/python3.12 '/mnt/c/Users/q287js/Documents/ChatGPT/Aspergillus genomics/gate_execution/2026-11-18_v1/scripts/resume_checkpoint.py' $resumeMode
if ($LASTEXITCODE -ne 0) { throw 'Checkpoint verification or resumption failed. Inspect output; do not force restart.' }
