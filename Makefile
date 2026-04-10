PYTHON ?= python
PIP ?= $(PYTHON) -m pip
STREAMLIT ?= $(PYTHON) -m streamlit
HOST ?= 0.0.0.0
PORT ?= 8501

.PHONY: setup pipeline dashboard

setup:
	$(PIP) install -r requirements.txt

pipeline:
	$(PYTHON) load_data.py
	$(PYTHON) generate_outputs.py

dashboard:
	$(STREAMLIT) run dashboard.py --server.address $(HOST) --server.port $(PORT)
