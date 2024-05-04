rm $1/*.txt
mv out_*.txt $1
for i in $1/*.txt
do
    ./transpose.sh $i > ${i}_tr.txt
done
tail -n 319 $1/out_C3.txt > $1/out_C3_cut.txt
tail -n +3 $1/out_C3.txt > $1/out_C3a.txt
rm $1/out_C3.txt
mv $1/out_C3a.txt $1/out_C3.txt
rm $1/out_C1.txt_tr.txt
rm $1/out_C2.txt_tr.txt
rm $1/out_C3.txt_tr.txt
rm $1/out_C0.txt_tr.txt

tail -n +13 $1/out_H.txt_tr.txt | head -n -2 > $1/a1.txt
rm $1/out_H.txt_tr.txt
mv $1/a1.txt $1/out_H.txt_tr.txt

tail -n +6 $1/out_T.txt_tr.txt > $1/a1.txt
rm $1/out_T.txt_tr.txt
mv $1/a1.txt $1/out_T.txt_tr.txt

tail -n +4 $1/out_V.txt_tr.txt > $1/a1.txt
rm $1/out_V.txt_tr.txt
mv $1/a1.txt $1/out_V.txt_tr.txt

tail -n +4 $1/out_ch4.txt_tr.txt > $1/a1.txt
rm $1/out_ch4.txt_tr.txt
mv $1/a1.txt $1/out_ch4.txt_tr.txt

tail -n +4 $1/out_co2.txt_tr.txt > $1/a1.txt
rm $1/out_co2.txt_tr.txt
mv $1/a1.txt $1/out_co2.txt_tr.txt

tail -n +4 $1/out_rho.txt_tr.txt > $1/a1.txt
rm $1/out_rho.txt_tr.txt
mv $1/a1.txt $1/out_rho.txt_tr.txt

tail -n +4 $1/out_soc.txt_tr.txt > $1/a1.txt
rm $1/out_soc.txt_tr.txt
mv $1/a1.txt $1/out_soc.txt_tr.txt

./yearly_avg.sh $1/out_C3.txt > $1/yearly_avg.txt