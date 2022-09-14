L(D)=
{
	my(d,l,u,f);
	[d,l]=coredisc(D,1);
	u=quadclassunit(d);
	f=factor(l);
	[D,d,l,(2*u.clgp[1]*u.reg/sqrt(D))*
	prod(i=1,#f[,1],1+(f[i,1]-kronecker(d,f[i,1]))*(f[i,1]^f[i,2]-1)/(f[i,1]-1))]
}

\\ run L(D) for all D <= x congruent to a mod q and output to file
\\ skips squares and assumes q is odd
go(a,q,x,outfile,prec)=
{
	my(fp,v);
	default(realprecision,prec);
	fp=fileopen(outfile,"w");
	for(amod4=0,1,
		forstep(D=lift(chinese(Mod(a,q),Mod(amod4,4))),x,4*q,
			if(!issquare(D),
				v=L(D);
				filewrite(fp,v[1]" "v[2]" "v[3]" "v[4]);
			);
		);
	);
}
