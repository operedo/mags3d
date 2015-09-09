#You must define COL with option -v of awk
# usage:
# cat file.txt | awk -v COL=1 -W exec=print_average.awk
{
if(min=="")
{
min=max=$COL
};
if($COL>max)
{
max=$COL
maxIdx=NR
};
if($COL< min)
{
min=$COL
minIdx=NR
};
total+=$COL;
count+=1
}
END {print total/count}
