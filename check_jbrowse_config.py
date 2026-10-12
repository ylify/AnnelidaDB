import os
import urllib.request
import urllib.parse
import urllib.error

root = 'static/data/Owenia_fusiformis/genome'
gffpath = None
for fname in os.listdir(root):
    if fname.endswith('.gff3') or fname.endswith('.gff3.gz'):
        gffpath = os.path.join(root, fname)
        break
print('gffpath', gffpath)
if not gffpath:
    raise SystemExit('no gff')
with open(gffpath, 'rt', encoding='utf-8', errors='ignore') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) > 4:
            seqid = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            loc = f'{seqid}:{start}..{end}'
            print('sample loc', loc)
            break
url = 'http://127.0.0.1:8001/jbrowse-config/?' + urllib.parse.urlencode({'data': '/static/data/Owenia_fusiformis/genome', 'loc': loc})
print('url', url)
try:
    with urllib.request.urlopen(url, timeout=10) as f:
        print(f.read().decode('utf-8')[:2000])
except urllib.error.HTTPError as e:
    print('HTTP ERROR', e.code)
    print(e.read().decode('utf-8', errors='replace'))
except Exception as e:
    print('ERROR', e)
