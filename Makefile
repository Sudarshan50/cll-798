.PHONY: data dev build docker-build docker-run test clean

# Full pipeline, in dependency order. 02 needs the 2 GB dataset at
# ../chemont_project/data/classyfire_main.tsv.zst and the `zstd` CLI; it takes
# ~7 minutes. 08 writes the 292 MB frontend/public/points.bin.
data:
	python3 scripts/01_verify_taxonomy.py
	python3 scripts/02_scan_dataset.py
	python3 scripts/03_verify_evaluation.py
	python3 scripts/05_build_graph_data.py
	python3 scripts/06_build_pointcloud.py
	python3 scripts/07_build_network.py
	python3 scripts/08_build_network_points.py

dev:
	cd frontend && npm install && npm run dev

build:
	cd frontend && npm ci && npm run build

docker-build:
	docker build -f deploy/Dockerfile -t chemical-space-atlas .

docker-run:
	docker compose -f deploy/docker-compose.yml up --build

test:
	python3 scripts/test_checks.py

clean:
	rm -rf frontend/dist frontend/public/points.bin frontend/public/layout.json
	find . -name __pycache__ -type d -prune -exec rm -rf {} +
