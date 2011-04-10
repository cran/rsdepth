rsmed <- function(pt, ...)
{
	if(!is.matrix(pt)) stop ("first argument pt must be a 2 dimensional matrix");
        median = .C("rs_med",
                    points=as.double(pt), 
                    md=double(2), 
                    sz=as.integer(length(pt)/2),
                    PACKAGE = "rsdepth")$md;
	return (median);

}

