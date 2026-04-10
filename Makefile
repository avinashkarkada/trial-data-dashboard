PYTHON ?= python
PIP ?= $(PYTHON) -m pip
STREAMLIT ?= $(PYTHON) -m streamlit

.PHONY: setup pipeline dashboard

setup:
	$(PIP) install --upgrade pip
	$(PIP) install -r requirements.txt

pipeline:
	$(PYTHON) load_data.py
	$(PYTHON) generate_outputs.py

dashboard:
	$(STREAMLIT) run dashboard.py --server.address 0.0.0.0 --server.port 8501
