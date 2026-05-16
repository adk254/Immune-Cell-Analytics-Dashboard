setup:
	python -m pip install -r requirements.txt

pipeline:
	python load_data.py
	python analyze.py

dashboard:
	python -m streamlit run dashboard.py