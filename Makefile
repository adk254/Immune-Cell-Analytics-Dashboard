setup:
	pip install -r requirements.txt

pipeline:
	python load_data.py
	python analyze.py

dashboard:
	streamlit run dashboard.py