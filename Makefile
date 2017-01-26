init:
	gunzip db/*.gz
	./rRNA_pipeline.py --test

python3:
	cp source_py3/* .
