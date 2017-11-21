init:
	gunzip db/*.gz
	./rRNA_pipeline.py --test
