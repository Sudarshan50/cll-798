.PHONY: data extension ui-data dev build docker-build docker-run test clean

PY ?= python3

# Full pipeline, in dependency order. 02 needs the 2 GB dataset at
# ../chemont_project/data/classyfire_main.tsv.zst and the `zstd` CLI; it takes
# ~7 minutes. 08 writes the 292 MB frontend/public/points.bin.
data:
	$(PY) scripts/01_verify_taxonomy.py
	$(PY) scripts/02_scan_dataset.py
	$(PY) scripts/03_verify_evaluation.py
	$(PY) scripts/05_build_graph_data.py
	$(PY) scripts/06_build_pointcloud.py
	$(PY) scripts/07_build_network.py
	$(PY) scripts/08_build_network_points.py

# Extension (Phase 2 + bipartite). 09 is the only step that reads the 2 GB dataset
# (~4 minutes, 6 workers); everything after it works off its sufficient statistics.
# Needs numpy, scipy, networkx, powerlaw, openpyxl: make extension PY=.venv/bin/python
extension:
	$(PY) scripts/09_scan_bipartite.py
	$(PY) scripts/09b_alt_cooccurrence.py
	$(PY) scripts/10_build_overlay.py
	$(PY) scripts/11_separate_edges.py
	$(PY) scripts/12_validate_flagged.py
	$(PY) scripts/13_network_measures.py
	$(PY) scripts/13b_error_propagation.py
	$(PY) scripts/13c_build_ui.py
	$(PY) scripts/14_bipartite_projections.py
	$(PY) scripts/15_bipartite_backbones.py
	$(PY) scripts/16_functional_groups.py
	$(MAKE) ui-data

ui-data:
	cp results/phase2/ui.json frontend/public/phase2.json
	cp results/bipartite/ui.json frontend/public/bipartite.json

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
	$(PY) scripts/test_phase2.py
	$(PY) scripts/test_bipartite.py

clean:
	rm -rf frontend/dist frontend/public/points.bin frontend/public/layout.json
	find . -name __pycache__ -type d -prune -exec rm -rf {} +
