# Deploying the Chemical Space Atlas

nginx serving a static Vite bundle. The two data assets are generated separately
and mounted in; they are not part of the image.

## 1. Generate the data assets

```bash
pip install -r requirements.txt
make data                     # full pipeline; needs the 2 GB dataset + `zstd`
```

If `results/` is already populated, only the last two steps are needed:

```bash
python3 scripts/07_build_network.py         # results/network_layout.json
python3 scripts/08_build_network_points.py  # frontend/public/{points.bin,layout.json}
```

This produces `frontend/public/points.bin` (292 MB) and
`frontend/public/layout.json` (~3.6 MB).

## 2. Build and run

```bash
docker compose -f deploy/docker-compose.yml up --build   # http://localhost:8080
```

Compose mounts `frontend/public/` read-only at `/srv/atlas-data`; nginx serves
`/points.bin` and `/layout.json` from there. To regenerate the data, rerun step 1
on the host and reload the page — no rebuild.

Without compose:

```bash
docker build -f deploy/Dockerfile -t chemical-space-atlas .
docker run -p 8080:80 -v "$PWD/frontend/public:/srv/atlas-data:ro" chemical-space-atlas
```

## Resources

**Server:** static files only. Any small instance is enough; the cost is 292 MB
of egress per full page load, so bandwidth, not CPU or RAM, is the limit.

**Client:** the browser must allocate a single 292 MB WebGL2 vertex buffer, so it
needs WebGL2 and a GPU with roughly 512 MB free. `src/gl.js` streams `points.bin`
and uploads it chunk by chunk, dropping each chunk, so peak CPU memory is one
chunk rather than 292 MB — but the GPU allocation is all-or-nothing.

**The failure mode:** on a machine with under about 8 GB of RAM (or an integrated
GPU sharing it), the buffer allocation or the stream fails partway. It is not
silent — the point counter in the corner stops climbing before 73,105,281, and
that number is the count actually resident on the GPU. If it stalls, the client
ran out of memory; there is nothing to fix server-side.
