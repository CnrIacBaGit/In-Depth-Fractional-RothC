awk '
BEGIN { 
    co2=0;
    ch4=0;
    co2_max=0;
    co2_min=1000;
    ch4_max=0;
    ch4_min=1000;
    year=0;
}
{
if (int($2/365.0)!=year)
{
    printf "%d %g %g %g %g %g %g\n",year,$127-co2,$129-ch4,co2_max,co2_min,ch4_max,ch4_min
    year=year+1;
    co2=$127;
    ch4=$129;
    co2_max=0;
    co2_min=1000;
    ch4_max=0;
    ch4_min=1000;
}
lco2=$127;
lch4=$129;
if ($117 > ch4_max) ch4_max=$117
if ($117 < ch4_min) ch4_min=$117
if ($113 > co2_max) co2_max=$113
if ($113 < co2_min) co2_min=$113
}
END {
    printf "%d %g %g %g %g %g %g\n",year,lco2-co2,lch4-ch4,co2_max,co2_min,ch4_max,ch4_min
}
' $1