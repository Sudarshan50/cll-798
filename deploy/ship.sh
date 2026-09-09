#!/usr/bin/env bash
# Ship the built atlas to a remote host. Only three things travel: the compiled app
# (~176 KB), layout.json, and points.bin. No sources, no node_modules, no dataset.
#
#   ./deploy/ship.sh ec2-user@host /path/to/key.pem
#
# points.bin is re-sent only when the remote checksum differs, so re-running this
# after a frontend-only change moves kilobytes, not 292 MB.
set -euo pipefail
HOST=${1:?usage: ship.sh user@host key.pem}
KEY=${2:?usage: ship.sh user@host key.pem}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
R() { ssh -i "$KEY" -o BatchMode=yes "$HOST" "$@"; }

[ -f "$ROOT/frontend/dist/index.html" ] || { echo "build first: cd frontend && npm run build"; exit 1; }

R 'mkdir -p ~/site/html ~/site/data'
tar -czf - -C "$ROOT/frontend/dist" index.html assets | R 'tar -xzf - -C ~/site/html'
tar -czf - -C "$ROOT/deploy" nginx.conf | R 'tar -xzf - -C ~/site'
R 'cat > ~/site/data/layout.json' < "$ROOT/frontend/public/layout.json"

want=$(shasum -a 256 "$ROOT/frontend/public/points.bin" | cut -d' ' -f1)
have=$(R 'sha256sum ~/site/data/points.bin 2>/dev/null | cut -d" " -f1' || true)
if [ "$want" != "$have" ]; then
  echo "shipping points.bin (292 MB)"
  R 'cat > ~/site/data/points.bin' < "$ROOT/frontend/public/points.bin"
  have=$(R 'sha256sum ~/site/data/points.bin | cut -d" " -f1')
  [ "$want" = "$have" ] || { echo "checksum mismatch after transfer"; exit 1; }
else
  echo "points.bin unchanged, skipped"
fi

R 'sudo docker rm -f web-2 >/dev/null 2>&1 || true
   sudo docker run -d --name web-2 --restart unless-stopped -p 80:80 \
     -v ~/site/html:/usr/share/nginx/html:ro \
     -v ~/site/data:/srv/atlas-data:ro \
     -v ~/site/nginx.conf:/etc/nginx/nginx.conf:ro \
     nginx:alpine >/dev/null && echo "web-2 running on :80"'
