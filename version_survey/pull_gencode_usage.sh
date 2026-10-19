#!/bin/bash
# GitHub commits mentioning each GENCODE version, as a proxy for adoption over time
set -euo pipefail

OUTDIR="$(dirname "$0")"
OUTFILE="$OUTDIR/gencode_version_commit_mentions.tsv"

VERSIONS=(25 30 35 39 44 45 46 47 48 49 50)
MAX_PAGES=3

echo -e "version\tdate\tsha\trepo" > "$OUTFILE"

for v in "${VERSIONS[@]}"; do
    query="gencode.v${v}"
    page=1
    while [ "$page" -le "$MAX_PAGES" ]; do
        resp=$(curl -s -G "https://api.github.com/search/commits" \
            --data-urlencode "q=${query}" \
            --data-urlencode "sort=author-date" \
            --data-urlencode "order=asc" \
            --data-urlencode "per_page=100" \
            --data-urlencode "page=${page}" \
            -H "Accept: application/vnd.github.cloak-preview+json")

        if echo "$resp" | python3 -c "import json,sys; d=json.load(sys.stdin); sys.exit(0 if 'items' in d else 1)" 2>/dev/null; then
            :
        else
            echo "ERROR on v${v} page ${page}: $resp" >&2
            sleep 15
            continue
        fi

        n=$(echo "$resp" | python3 -c "import json,sys; d=json.load(sys.stdin); print(len(d.get('items', [])))")
        if [ "$n" = "0" ]; then
            break
        fi

        echo "$resp" | python3 -c "
import json, sys
d = json.load(sys.stdin)
for c in d.get('items', []):
    date = c['commit']['author']['date']
    sha = c['sha'][:8]
    repo = c['repository']['full_name']
    print(f'${v}\t{date}\t{sha}\t{repo}')
" >> "$OUTFILE"

        page=$((page + 1))
        sleep 12  # 10 requests/min limit
    done
    echo "Pulled v${v}: done" >&2
done

echo "Done. Output: $OUTFILE" >&2
