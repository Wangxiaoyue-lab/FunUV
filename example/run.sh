python3 ../Sub_MutSeqByMANE.py -v test_5utr -r Homo_sapiens_assembly38.fasta > test_5utr.addseq.txt
python3 ../Sub_MutSeqByMANE.py -v test_3utr -r Homo_sapiens_assembly38.fasta > test_3utr.addseq.txt

python3 ../FunUV_annotation.py -f test_5utr.addseq.txt -o test_5utr -r 5UTR
python3 ../FunUV_annotation.py -f test_3utr.addseq.txt -o test_3utr -r 3UTR

python3 ../Prepare_forPred.py -v test_5utr_annotation.txt -r 5UTR > test_5utr_subanno.txt
python3 ../Prepare_forPred.py -v test_3utr_annotation.txt -r 3UTR > test_3utr_subanno.txt

python3 ../FunUV_predict.py -pf test_5utr_subanno.txt -r 5UTR -o test_5utr
python3 ../FunUV_predict.py -pf test_3utr_subanno.txt -r 3UTR -o test_3utr
