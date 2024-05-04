awk '
BEGIN { 
    co2=0;
    ch4=0;
    year=0;
}
{
if (int($2/365.0)!=year)
{
    printf "%d %g %g\n",year,$127-co2,$129-ch4
    year=year+1;
    co2=$127;
    ch4=$129;
}
lco2=$127;
lch4=$129;
}
END {
    printf "%d %g %g\n",year,lco2-co2,lch4-ch4
}
' $1