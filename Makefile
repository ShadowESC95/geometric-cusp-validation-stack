.PHONY: install fast full test examples amp gallery manifest backends clean
install:
	python -m pip install -e '.[all,test]'
fast:
	PYTHONPATH=src python scripts/run_full_research.py --fast
full:
	PYTHONPATH=src python scripts/run_full_research.py
amp:
	PYTHONPATH=src python scripts/run_amplituhedron_sweep.py

gallery:
	PYTHONPATH=src python scripts/make_gallery.py
manifest:
	PYTHONPATH=src python scripts/build_manifest.py
backends:
	PYTHONPATH=src python scripts/backend_report.py

test:
	PYTHONPATH=src pytest -q
examples:
	PYTHONPATH=src python -m examples.run_suite
clean:
	find . -type d -name '__pycache__' -prune -exec rm -rf {} +
