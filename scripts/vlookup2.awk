BEGIN{OFS="\t"}
FNR==NR{
  a[$1]=$1
  next
}
{  if ($12 in a) {print $0, "RM"} else {print $0, "KP"} }
