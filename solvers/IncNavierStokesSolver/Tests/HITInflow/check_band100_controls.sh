#!/usr/bin/env bash
# Read-only qualification progress from the Mac or Linux; no solver launch.
set -euo pipefail
hit_control_root=${1:?Usage: check_band100_controls.sh LINUX_CONTROL_QUEUE_DIRECTORY}
if [[ $(uname -s) == Linux ]]; then
    hit_status_command=(python3 - "$hit_control_root")
else
    printf -v hit_control_quoted '%q' "$hit_control_root"
    hit_status_command=(ssh -o BatchMode=yes -o ConnectTimeout=15 jay@100.105.247.110 "python3 - $hit_control_quoted")
fi
"${hit_status_command[@]}" <<'PY'
import datetime,json,re,sys
from pathlib import Path
root=Path(sys.argv[1]);queue=json.loads((root/'queue.json').read_text())
print(datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')+' — CBC100 qualification')
print('Queue: '+queue['status'])
hold_path=root.parent/'resource-hold-01.json'
if hold_path.is_file():
    hold=json.loads(hold_path.read_text())
    if hold.get('queue_pid')==queue.get('pid'):
        print('Resource scheduling: '+hold['status'])
        if hold['status']=='WAITING_FOR_EXISTING_GUARDS':
            print('HIT processes are held; automatic resume follows completion of existing server guards.')
        if hold.get('updated_utc'):print('Resource snapshot: '+hold['updated_utc'])
print('Completed controls: '+(', '.join(queue.get('completed_cases',[])) or 'none'))
current=queue.get('current')
if current:
    print('Stage: '+current['stage'])
    path=Path(current['output'])/'progress.json'
    if path.is_file():
        progress=json.loads(path.read_text());native=progress.get('current')
        if native:
            print('Command: '+native['command'])
            log=Path(native['log'])
            if log.is_file():
                with log.open('rb') as stream:
                    stream.seek(max(0,log.stat().st_size-131072));text=stream.read().decode(errors='replace')
                rows=re.findall(r'Steps:\s*(\d+)\s+Time:\s*([\d.eE+-]+)',text)
                if rows:
                    step,t=rows[-1];total=native.get('total_steps')
                    print(f'Step: {step} / {total or "?"}; simulated time: {float(t):.6f}s')
                    if total:print(f'Progress: {100*int(step)/total:.2f}%')
                else:print('No step report yet: native initialization or first info interval.')
        if progress.get('updated_utc'):print('Solver snapshot: '+progress['updated_utc'])
        if 'budget_gate' in progress:
            gate=progress['budget_gate'];print(f"Budget gate: {'pass' if gate['passed'] else 'reject'}; forecast {gate['forecast_wall_seconds']/3600:.2f}h, sampled RSS {gate['sampled_peak_tree_rss_gib']:.2f}GiB")
if queue.get('error'):print('Error: '+queue['error'])
print('Queue stage updated: '+queue['updated_utc'])
print('Full-distance production: HELD pending qualification review')
PY
