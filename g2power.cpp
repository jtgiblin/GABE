#include "g2header.h"

#define lmINDEX(l,m) l*(l+1)+m

void deIndex(gIdx num, gIdx *i, gIdx *j, gIdx *k)
{
	*i=num/N/N;
	*j=(num%(N*N))/N;
	*k=(num%(N*N))%N;
}

gIdx countSphere(gIdx radSphere){
	gIdx count=0;
	gIdx i,j,k;
	
	LOOP{
		if(
			(i-N/2)*(i-N/2)+(k-N/2)*(k-N/2)+(j-N/2)*(j-N/2)<(radSphere-1)*(radSphere-1)
			&&(i-N/2)*(i-N/2)+(k-N/2)*(k-N/2)+(j-N/2)*(j-N/2)>=(radSphere-2)*(radSphere-2)
			){
			count++;
		}
	}
	return count;
}

void makeSphere(gNum sphere[], gIdx nSphere, gIdx radSphere){
	gIdx count=0;
	gIdx i,j,k;
	
	LOOP{
		if(
			(i-N/2)*(i-N/2)+(k-N/2)*(k-N/2)+(j-N/2)*(j-N/2)<(radSphere-1)*(radSphere-1)
			&&(i-N/2)*(i-N/2)+(k-N/2)*(k-N/2)+(j-N/2)*(j-N/2)>=(radSphere-2)*(radSphere-2)
			){
			sphere[count]=k+N*j+N*N*i;
			count++;
		}
	}
}

gNum rSphY(gIdx lmIdx, gNum cp, gNum sp, gNum ct, gNum st){
	switch(lmIdx)
	{
		//l=0
		case 0: 
			return 1./sqrtl(M_PI)/2.;
			break;
		//l=1
		case 1:
			return sqrtl(3./2./M_PI)*cp*st/2.;
			break;
		case 2:
			return sqrtl(3./M_PI)*ct/2.;
			break;
		case 3:
			return -sqrtl(3./2./M_PI)*cp*st/2.;
			break;
		//l=2
		case 4:
			return -(sqrtl(15./(2.*M_PI))*(cp - sp)*(cp + sp)*(-1. + ct*ct - st*st))/8.;
			break;
		case 5:
			return (cp*ct*sqrtl(15./(2.*M_PI))*st)/2.;
			break;
		case 6:
			return ((-1. + 3.*ct*ct)*sqrtl(5./M_PI))/4.;
			break;
		case 7:
			return -(cp*ct*sqrtl(15./(2.*M_PI))*st)/2.;
			break;
		case 8:
			return -(sqrtl(15./(2.*M_PI))*(cp - sp)*(cp + sp)*(-1. + ct*ct - st*st))/8.;
			break;
		//l=3
		case 9:
			return (cp*sqrtl(35./M_PI)*(cp*cp - 3.*sp*sp)*st*(3. - 3.*ct*ct + st*st))/32.;
			break;
		case 10:
			return -(ct*sqrtl(105./(2.*M_PI))*(cp - sp)*(cp + sp)*(-1. + ct*ct - 3.*st*st))/16.;
			break;
		case 11:
			return (cp*sqrtl(21./M_PI)*st*(1. + 15.*ct*ct - 5.*st*st))/32.;
			break;
		case 12:
			return (ct*sqrtl(7./M_PI)*(3. + 5.*ct*ct - 15.*st*st))/16.;
			break;
		case 13:
			return (cp*sqrtl(21./M_PI)*st*(-1. - 15.*ct*ct + 5.*st*st))/32.;
			break;
		case 14:
			return -(ct*sqrtl(105./(2.*M_PI))*(cp - sp)*(cp + sp)*(-1. + ct*ct - 3.*st*st))/16.;
			break;
		case 15:
			return -(cp*sqrtl(35./M_PI)*(cp*cp - 3.*sp*sp)*st*(3. - 3.*ct*ct + st*st))/32.;
			break;
	}
}
gNum cSphY(gIdx lmIdx, gNum cp, gNum sp, gNum ct, gNum st){
	switch(lmIdx)
	{
		//l=0
		case 0:
			return 0;
			break;
		//l=1
		case 1:
			return -(sqrtl(3./(2.*M_PI))*sp*st)/2.;
			break;
		case 2:
			return 0;
			break;
		case 3:
			return -(sqrtl(3./(2.*M_PI))*sp*st)/2.;
			break;
		//l=2
		case 4:
			return (cp*sqrtl(15./(2.*M_PI))*sp*(-1. + ct*ct - st*st))/4.;
			break;
		case 5:
			return -(ct*sqrtl(15./(2.*M_PI))*sp*st)/2.;
			break;
		case 6:
			return 0;
			break;
		case 7:
			return -(ct*sqrtl(15./(2.*M_PI))*sp*st)/2.;
			break;
		case 8:
			return (cp*sqrtl(15./(2.*M_PI))*sp*(1. - ct*ct + st*st))/4.;
			break;
		//l=3
		case 9:
			return (sqrtl(35./M_PI)*sp*(-3.*cp*cp + sp*sp)*st*(3. - 3.*ct*ct + st*st))/32.;
			break;
		case 10:
			return (cp*ct*sqrtl(105./(2.*M_PI))*sp*(-1. + ct*ct - 3.*st*st))/8.;
			break;
		case 11:
			return (sqrtl(21./M_PI)*sp*st*(-1. - 15.*ct*ct + 5.*st*st))/32.;
			break;
		case 12:
			return 0;
			break;
		case 13:
			return (sqrtl(21./M_PI)*sp*st*(-1. - 15.*ct*ct + 5.*st*st))/32.;
			break;
		case 14:
			return -(cp*ct*sqrtl(105./(2.*M_PI))*sp*(-1. + ct*ct - 3.*st*st))/8.;
			break;
		case 15:
			return (sqrtl(35./M_PI)*sp*(-3.*cp*cp + sp*sp)*st*(3. - 3.*ct*ct + st*st))/32.;
			break;
	}
}

void modePowerOut(gNum tin, int first, gIdx RR){
	gIdx idx,i,j,k;
	gNum x,y,z,rPi,tPi;
	gNum cp,sp,ct,st,totPow;
	//the size of these are 16 to cary upto octopole terms
	gNum almc[16]={0.};
	gNum almr[16]={0.};
	gNum blmc[16]={0.};
	gNum blmr[16]={0.};
	gIdx nSphere=countSphere(RR);
	gNum sphere[nSphere];//this stores the index of the points along our sphere
	makeSphere(sphere,nSphere,RR);
	const gNum dOmega=4.*M_PI/(double)nSphere;
	totPow=0.;
	for(idx=0;idx<nSphere;idx++){
		deIndex(sphere[idx],&i,&j,&k);
		rPi=dfdr(field,0,0,i,j,k)- dfdr_analytic(0,0,i,j,k,tin);
		tPi=dfield[INDEX(0,0,i,j,k)];
		x=((double)i-(double)N/2.);
		y=((double)j-(double)N/2.);
		z=((double)k-(double)N/2.);
		if(sqrtl((x*x)+(y*y))==0){
			cp=0.;
			sp=0.;
		}
		else{
			cp=y/sqrtl((x*x)+(y*y));
			sp=x/sqrtl((x*x)+(y*y));
		}
		ct=z/sqrtl((x*x)+(y*y)+(z*z));
		st=sqrtl((x*x)+(y*y))/sqrtl((x*x)+(y*y)+(z*z));
		for(gIdx l=0;l<4;l++){
			for(gIdx m=-l;m<=l;m++) {
				almc[lmINDEX(l,m)]+=st*rPi*cSphY(lmINDEX(l,m),cp,sp,ct,st);
				almr[lmINDEX(l,m)]+=st*rPi*rSphY(lmINDEX(l,m),cp,sp,ct,st);
				blmc[lmINDEX(l,m)]+=st*tPi*cSphY(lmINDEX(l,m),cp,sp,ct,st);
				blmr[lmINDEX(l,m)]+=st*tPi*rSphY(lmINDEX(l,m),cp,sp,ct,st);
			}
		}
		totPow+=st*rPi*tPi;
	}

	static FILE *powerout;
	char powername[100];
	sprintf(powername, "./slices/power%d.dat",RR);
	powerout=fopen(powername,"a");
	if((totPow!=0. && totPow/totPow!=1.)) // These two separate checks between them work on all the compilers I've tested
    {
        printf("Unstable solution developed. Pi not numerical at t=%Le\n",tin);
        output_parameters();
        exit(1);
    }

    fprintf(powerout, "%Le ",totPow*dOmega*(RR-.5)*(RR-.5)*dx*dx);
	for(gIdx l=0;l<4;l++){
		for(gIdx m=-l;m<=l;m++) {
		fprintf(powerout, "%Le %Le %Le %Le ",almc[lmINDEX(l,m)]*dOmega,almr[lmINDEX(l,m)]*dOmega,blmc[lmINDEX(l,m)]*dOmega,blmr[lmINDEX(l,m)]*dOmega);
		}
	}
	fprintf(powerout, "\n");
	fclose(powerout);
}