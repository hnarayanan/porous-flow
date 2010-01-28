function truncate(x)
{
    EPS=1.e-16
    absx = x >= 0 ? x : -x
    if (absx < EPS) {
	return "0.0"
    }
    else {
	return x
    }
}

{
    if (NF == 1)
	$i = truncate($i);
    print
}