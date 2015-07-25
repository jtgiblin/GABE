#include "g2header.h"
#include "g2sphericalharmonics.h"

#define lmINDEX(l,m) l*(l+1)+m
#define lmCFINDEX(l,m) l+m

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

void modePowerOut(gNum tin, int first, gIdx RR){
	gIdx idx,i,j,k;
	gNum x,y,z,rPi,tPi;
	gNum cp,sp,ct,st,totPow;
	//the size of these are 16 to cary upto octopole terms
	gNum almc[(MAX_MODE+1)*(MAX_MODE+1)]={0.};
	gNum almr[(MAX_MODE+1)*(MAX_MODE+1)]={0.};
	gNum blmc[(MAX_MODE+1)*(MAX_MODE+1)]={0.};
	gNum blmr[(MAX_MODE+1)*(MAX_MODE+1)]={0.};
	gNum plr[MAX_MODE+1]={0.};
	gNum plc[MAX_MODE+1]={0.};
	gNum p[2]={0.};
	gIdx nSphere=countSphere(RR);
	gNum sphere[nSphere];//this stores the index of the points along our sphere
	makeSphere(sphere,nSphere,RR);
	gNum dOmega=4.*M_PI/(double)nSphere;
	gNum rSquared=(RR-.5)*(RR-.5)*dx*dx;
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
		for(gIdx l=0;l<=MAX_MODE;l++){
			for(gIdx m=-l;m<=l;m++) {
				almc[lmINDEX(l,m)]+=rPi*cSphY(lmINDEX(l,m),cp,sp,ct,st);
				almr[lmINDEX(l,m)]+=rPi*rSphY(lmINDEX(l,m),cp,sp,ct,st);
				blmc[lmINDEX(l,m)]+=tPi*cSphY(lmINDEX(l,m),cp,sp,ct,st);
				blmr[lmINDEX(l,m)]+=tPi*rSphY(lmINDEX(l,m),cp,sp,ct,st);
			}
		}
		totPow+=rPi*tPi;
	}

	for(gIdx l=0;l<=MAX_MODE;l++){
		for(gIdx m=-l;m<=l;m++) {
			plr[l]+=almc[lmINDEX(l,m)]*blmc[lmINDEX(l,m)]+almr[lmINDEX(l,m)]*blmr[lmINDEX(l,m)];
			plc[l]+=almr[lmINDEX(l,m)]*blmc[lmINDEX(l,m)]-almc[lmINDEX(l,m)]*blmr[lmINDEX(l,m)];
		}
		p[0]+=plr[l];
		p[1]+=plc[l];
	}


	static FILE *powerout;
	char powername[100];
	sprintf(powername, "./slices/power%d.dat",RR);
	powerout=fopen(powername,"a");
	if((totPow!=0. && totPow/totPow!=1.)) // These two separate checks between them work on all the compilers we have used
    {
        printf("Unstable solution developed. Pi not numerical at t=%Le\n",tin);
        output_parameters();
        exit(1);
    }

    fprintf(powerout, "%Le %Le %Le ",totPow*dOmega*rSquared,p[0]*dOmega*dOmega*rSquared,p[1]*dOmega*dOmega*rSquared);
	for(gIdx l=0;l<=MAX_MODE;l++){
		fprintf(powerout, "%Le %Le ",plr[l]*dOmega*dOmega*rSquared,plc[l]*dOmega*dOmega*rSquared);
		
		//For use if you want all the coefficients
		/*
		for(gIdx m=-l;m<=l;m++) {
			fprintf(powerout, "%Le %Le %Le %Le ",almc[lmINDEX(l,m)]*dOmega,almr[lmINDEX(l,m)]*dOmega,blmc[lmINDEX(l,m)]*dOmega,blmr[lmINDEX(l,m)]*dOmega);
		}
		*/
	}
	fprintf(powerout, "\n");
	fclose(powerout);
}