FILE=opt_12_g3

source ../../groups.sh

for f in ${group3[@]}
do
    echo -n $f" " >> $FILE
done
echo "" >> $FILE
echo "opt"
cat opt.out | grep Time | awk 'BEGIN {ORS=" "; sum=0} {sum=sum+$3} (NR%10==0) {print sum/10; sum=0} END {print "\n"}' >> $FILE
echo "D100"
cat D100.out | grep Time | awk 'BEGIN {ORS=" "; sum=0} {sum=sum+$3} (NR%10==0) {print sum/10; sum=0} END {print "\n"}' >> $FILE
echo "D300"
cat D300.out | grep Time | awk 'BEGIN {ORS=" "; sum=0} {sum=sum+$3} (NR%10==0) {print sum/10; sum=0} END {print "\n"}' >> $FILE
echo "D500"
cat D500.out | grep Time | awk 'BEGIN {ORS=" "; sum=0} {sum=sum+$3} (NR%10==0) {print sum/10; sum=0} END {print "\n"}' >> $FILE
echo "D700"
cat D700.out | grep Time | awk 'BEGIN {ORS=" "; sum=0} {sum=sum+$3} (NR%10==0) {print sum/10; sum=0} END {print "\n"}' >> $FILE
echo "S50k"
cat S50k.out | grep Time | awk 'BEGIN {ORS=" "; sum=0} {sum=sum+$3} (NR%10==0) {print sum/10; sum=0} END {print "\n"}' >> $FILE
echo "S70k"
cat S70k.out | grep Time | awk 'BEGIN {ORS=" "; sum=0} {sum=sum+$3} (NR%10==0) {print sum/10; sum=0} END {print "\n"}' >> $FILE
echo "S90k"
cat S90k.out | grep Time | awk 'BEGIN {ORS=" "; sum=0} {sum=sum+$3} (NR%10==0) {print sum/10; sum=0} END {print "\n"}' >> $FILE
echo "nonested"
cat nonested.out | grep Time | awk 'BEGIN {ORS=" "; sum=0} {sum=sum+$3} (NR%10==0) {print sum/10; sum=0} END {print "\n"}' >> $FILE
echo "nested"
cat nested.out | grep Time | awk 'BEGIN {ORS=" "; sum=0} {sum=sum+$3} (NR%10==0) {print sum/10; sum=0} END {print "\n"}' >> $FILE

for f in ${group3[@]}
do
    echo -n $f" " >> $FILE"_MAX"
done
echo "" >> $FILE"_MAX"
echo "opt"
cat opt.out | grep Time | awk 'BEGIN {ORS=" "; sum=0} {sum=(sum>$3?sum:$3)} (NR%10==0) {print sum; sum=0} END {print "\n"}' >> $FILE"_MAX"
echo "D100"
cat D100.out | grep Time | awk 'BEGIN {ORS=" "; sum=0} {sum=(sum>$3?sum:$3)} (NR%10==0) {print sum; sum=0} END {print "\n"}' >> $FILE"_MAX"
echo "D300"
cat D300.out | grep Time | awk 'BEGIN {ORS=" "; sum=0} {sum=(sum>$3?sum:$3)} (NR%10==0) {print sum; sum=0} END {print "\n"}' >> $FILE"_MAX"
echo "D500"
cat D500.out | grep Time | awk 'BEGIN {ORS=" "; sum=0} {sum=(sum>$3?sum:$3)} (NR%10==0) {print sum; sum=0} END {print "\n"}' >> $FILE"_MAX"
echo "D700"
cat D700.out | grep Time | awk 'BEGIN {ORS=" "; sum=0} {sum=(sum>$3?sum:$3)} (NR%10==0) {print sum; sum=0} END {print "\n"}' >> $FILE"_MAX"
echo "S50k"
cat S50k.out | grep Time | awk 'BEGIN {ORS=" "; sum=0} {sum=(sum>$3?sum:$3)} (NR%10==0) {print sum; sum=0} END {print "\n"}' >> $FILE"_MAX"
echo "S70k"
cat S70k.out | grep Time | awk 'BEGIN {ORS=" "; sum=0} {sum=(sum>$3?sum:$3)} (NR%10==0) {print sum; sum=0} END {print "\n"}' >> $FILE"_MAX"
echo "S90k"
cat S90k.out | grep Time | awk 'BEGIN {ORS=" "; sum=0} {sum=(sum>$3?sum:$3)} (NR%10==0) {print sum; sum=0} END {print "\n"}' >> $FILE"_MAX"
echo "nonested"
cat nonested.out | grep Time | awk 'BEGIN {ORS=" "; sum=0} {sum=(sum>$3?sum:$3)} (NR%10==0) {print sum; sum=0} END {print "\n"}' >> $FILE"_MAX"
echo "nested"
cat nested.out | grep Time | awk 'BEGIN {ORS=" "; sum=0} {sum=(sum>$3?sum:$3)} (NR%10==0) {print sum; sum=0} END {print "\n"}' >> $FILE"_MAX"

for f in ${group3[@]}
do
    echo -n $f" " >> $FILE"_MIN"
done
echo "" >> $FILE"_MIN"
echo "opt"
cat opt.out | grep Time | awk 'BEGIN {ORS=" "; sum=100000} {sum=(sum>$3?$3:sum)} (NR%10==0) {print sum; sum=100000} END {print "\n"}' >> $FILE"_MIN"
echo "D100"
cat D100.out | grep Time | awk 'BEGIN {ORS=" "; sum=100000} {sum=(sum>$3?$3:sum)} (NR%10==0) {print sum; sum=100000} END {print "\n"}' >> $FILE"_MIN"
echo "D300"
cat D300.out | grep Time | awk 'BEGIN {ORS=" "; sum=100000} {sum=(sum>$3?$3:sum)} (NR%10==0) {print sum; sum=100000} END {print "\n"}' >> $FILE"_MIN"
echo "D500"
cat D500.out | grep Time | awk 'BEGIN {ORS=" "; sum=100000} {sum=(sum>$3?$3:sum)} (NR%10==0) {print sum; sum=100000} END {print "\n"}' >> $FILE"_MIN"
echo "D700"
cat D700.out | grep Time | awk 'BEGIN {ORS=" "; sum=100000} {sum=(sum>$3?$3:sum)} (NR%10==0) {print sum; sum=100000} END {print "\n"}' >> $FILE"_MIN"
echo "S50k"
cat S50k.out | grep Time | awk 'BEGIN {ORS=" "; sum=100000} {sum=(sum>$3?$3:sum)} (NR%10==0) {print sum; sum=100000} END {print "\n"}' >> $FILE"_MIN"
echo "S70k"
cat S70k.out | grep Time | awk 'BEGIN {ORS=" "; sum=100000} {sum=(sum>$3?$3:sum)} (NR%10==0) {print sum; sum=100000} END {print "\n"}' >> $FILE"_MIN"
echo "S90k"
cat S90k.out | grep Time | awk 'BEGIN {ORS=" "; sum=100000} {sum=(sum>$3?$3:sum)} (NR%10==0) {print sum; sum=100000} END {print "\n"}' >> $FILE"_MIN"
echo "nonested"
cat nonested.out | grep Time | awk 'BEGIN {ORS=" "; sum=100000} {sum=(sum>$3?$3:sum)} (NR%10==0) {print sum; sum=100000} END {print "\n"}' >> $FILE"_MIN"
echo "nested"
cat nested.out | grep Time | awk 'BEGIN {ORS=" "; sum=100000} {sum=(sum>$3?$3:sum)} (NR%10==0) {print sum; sum=100000} END {print "\n"}' >> $FILE"_MIN"

for f in ${group3[@]}
do
    echo -n $f" " >> $FILE"_T"
done
echo "" >> $FILE"_T"
echo "opt"
cat opt.out | grep tasks | awk 'BEGIN {ORS=" "; sum=0} {sum=sum+$2} (NR%10==0) {print sum/10; sum=0} END {print "\n"}' >> $FILE"_T"
echo "D100"
cat D100.out | grep tasks | awk 'BEGIN {ORS=" "; sum=0} {sum=sum+$2} (NR%10==0) {print sum/10; sum=0} END {print "\n"}' >> $FILE"_T"
echo "D300"
cat D300.out | grep tasks | awk 'BEGIN {ORS=" "; sum=0} {sum=sum+$2} (NR%10==0) {print sum/10; sum=0} END {print "\n"}' >> $FILE"_T"
echo "D500"
cat D500.out | grep tasks | awk 'BEGIN {ORS=" "; sum=0} {sum=sum+$2} (NR%10==0) {print sum/10; sum=0} END {print "\n"}' >> $FILE"_T"
echo "D700"
cat D700.out | grep tasks | awk 'BEGIN {ORS=" "; sum=0} {sum=sum+$2} (NR%10==0) {print sum/10; sum=0} END {print "\n"}' >> $FILE"_T"
echo "S50k"
cat S50k.out | grep tasks | awk 'BEGIN {ORS=" "; sum=0} {sum=sum+$2} (NR%10==0) {print sum/10; sum=0} END {print "\n"}' >> $FILE"_T"
echo "S70k"
cat S70k.out | grep tasks | awk 'BEGIN {ORS=" "; sum=0} {sum=sum+$2} (NR%10==0) {print sum/10; sum=0} END {print "\n"}' >> $FILE"_T"
echo "S90k"
cat S90k.out | grep tasks | awk 'BEGIN {ORS=" "; sum=0} {sum=sum+$2} (NR%10==0) {print sum/10; sum=0} END {print "\n"}' >> $FILE"_T"
echo "nonested"
cat nonested.out | grep tasks | awk 'BEGIN {ORS=" "; sum=0} {sum=sum+$2} (NR%10==0) {print sum/10; sum=0} END {print "\n"}' >> $FILE"_T"
echo "nested"
cat nested.out | grep tasks | awk 'BEGIN {ORS=" "; sum=0} {sum=sum+$2} (NR%10==0) {print sum/10; sum=0} END {print "\n"}' >> $FILE"_T"
