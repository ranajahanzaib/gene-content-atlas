PYTHON ?= python3
HOST ?= 127.0.0.1
PORT ?= 8000
export PANAROO_TIMEOUT_SECONDS ?= 86400

.PHONY: setup dev doctor seed-demo ingest-real analyze test lint build screenshots
setup:
	$(PYTHON) --version

dev:
	@HOST=$(HOST) PORT=$(PORT) $(PYTHON) services/api/app.py

doctor:
	@$(PYTHON) scripts/doctor.py

seed-demo:
	@$(PYTHON) scripts/seed_demo.py

ingest-real:
	@$(PYTHON) scripts/ingest_real.py

analyze:
	@$(PYTHON) scripts/run_panaroo.py

test:
	PYTHONPATH=services/api $(PYTHON) -m unittest discover -s services/api/tests -v

lint:
	$(PYTHON) -m compileall -q services/api scripts
	@$(PYTHON) scripts/check_frontend.py

build: lint
	@echo "Build checks passed. Static web assets are ready."

screenshots:
	@$(PYTHON) scripts/capture_screenshots.py
