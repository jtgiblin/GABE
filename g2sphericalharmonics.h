gNum rSphY(gIdx lmIdx, gNum cp, gNum sp, gNum ct, gNum st){
	gNum SP[7]={1.};
	gNum ST[7]={1.};
	gNum CP[7]={1.};
	gNum CT[7]={1.};
	for(int i=1;i<=6;i++){
		SP[i]=sp*SP[i-1];
		ST[i]=st*ST[i-1];
		CP[i]=cp*CP[i-1];
		CT[i]=ct*CT[i-1];
	}

	switch(lmIdx)
	{
	//l=0
		case 0:
			return 1/(2.*sqrtl(M_PI));
			break;
	
	//l=1
		case 1:
			return (cp*sqrtl(3/(2.*M_PI))*st)/2.;
			break;
		case 2:
			return (ct*sqrtl(3/M_PI))/2.;
			break;
		case 3:
			return -(cp*sqrtl(3/(2.*M_PI))*st)/2.;
			break;
	
	//l=2
		case 4:
			return -(sqrtl(15/(2.*M_PI))*(cp - sp)*(cp + sp)*(-1 + CT[2] - ST[2]))/8.;
			break;
		case 5:
			return (cp*ct*sqrtl(15/(2.*M_PI))*st)/2.;
			break;
		case 6:
			return ((-1 + 3*CT[2])*sqrtl(5/M_PI))/4.;
			break;
		case 7:
			return -(cp*ct*sqrtl(15/(2.*M_PI))*st)/2.;
			break;
		case 8:
			return -(sqrtl(15/(2.*M_PI))*(cp - sp)*(cp + sp)*(-1 + CT[2] - ST[2]))/8.;
			break;
	
	//l=3
		case 9:
			return (cp*sqrtl(35/M_PI)*(CP[2] - 3*SP[2])*st*(3 - 3*CT[2] + ST[2]))/32.;
			break;
		case 10:
			return -(ct*sqrtl(105/(2.*M_PI))*(cp - sp)*(cp + sp)*(-1 + CT[2] - 3*ST[2]))/16.;
			break;
		case 11:
			return (cp*sqrtl(21/M_PI)*st*(1 + 15*CT[2] - 5*ST[2]))/32.;
			break;
		case 12:
			return (ct*sqrtl(7/M_PI)*(3 + 5*CT[2] - 15*ST[2]))/16.;
			break;
		case 13:
			return (cp*sqrtl(21/M_PI)*st*(-1 - 15*CT[2] + 5*ST[2]))/32.;
			break;
		case 14:
			return -(ct*sqrtl(105/(2.*M_PI))*(cp - sp)*(cp + sp)*(-1 + CT[2] - 3*ST[2]))/16.;
			break;
		case 15:
			return -(cp*sqrtl(35/M_PI)*(CP[2] - 3*SP[2])*st*(3 - 3*CT[2] + ST[2]))/32.;
			break;
	
	//l=4
		case 16:
			return (3*sqrtl(35/(2.*M_PI))*(CP[4] - 6*CP[2]*SP[2] + SP[4])*(3 + CT[4] + 4*ST[2] + ST[4] - 2*CT[2]*(2 + 3*ST[2])))/128.;
			break;
		case 17:
			return (-3*cp*ct*sqrtl(35/M_PI)*(CP[2] - 3*SP[2])*st*(-1 + CT[2] - ST[2]))/16.;
			break;
		case 18:
			return (-3*sqrtl(5/(2.*M_PI))*(cp - sp)*(cp + sp)*(-3 + 7*CT[4] + 4*ST[2] + 7*ST[4] - 2*CT[2]*(2 + 21*ST[2])))/64.;
			break;
		case 19:
			return (3*cp*ct*sqrtl(5/M_PI)*st*(1 + 7*CT[2] - 7*ST[2]))/16.;
			break;
		case 20:
			return (3*(9 + 20*CT[2] + 35*CT[4] - 10*(2 + 21*CT[2])*ST[2] + 35*ST[4]))/(128.*sqrtl(M_PI));
			break;
		case 21:
			return (-3*cp*ct*sqrtl(5/M_PI)*st*(1 + 7*CT[2] - 7*ST[2]))/16.;
			break;
		case 22:
			return (-3*sqrtl(5/(2.*M_PI))*(cp - sp)*(cp + sp)*(-3 + 7*CT[4] + 4*ST[2] + 7*ST[4] - 2*CT[2]*(2 + 21*ST[2])))/64.;
			break;
		case 23:
			return (3*cp*ct*sqrtl(35/M_PI)*(CP[2] - 3*SP[2])*st*(-1 + CT[2] - ST[2]))/16.;
			break;
		case 24:
			return (3*sqrtl(35/(2.*M_PI))*(CP[4] - 6*CP[2]*SP[2] + SP[4])*(3 + CT[4] + 4*ST[2] + ST[4] - 2*CT[2]*(2 + 3*ST[2])))/128.;
			break;

	//l=5	
		case 25:
			return (3*cp*sqrtl(77/M_PI)*(CP[4] - 10*CP[2]*SP[2] + 5*SP[4])*st*(10 + 5*CT[4] + 5*ST[2] + ST[4] - 5*CT[2]*(3 + 2*ST[2])))/512.;
			break;
		case 26:
			return (3*ct*sqrtl(385/(2.*M_PI))*(CP[4] - 6*CP[2]*SP[2] + SP[4])*(2 + CT[4] + 9*ST[2] + 5*ST[4] - CT[2]*(3 + 10*ST[2])))/256.;
			break;
		case 27:
			return -(cp*sqrtl(385/M_PI)*(CP[2] - 3*SP[2])*st*(-6 + 45*CT[4] + 13*ST[2] + 9*ST[4] - 3*CT[2]*(13 + 30*ST[2])))/512.;
			break;
		case 28:
			return -(ct*sqrtl(1155/(2.*M_PI))*(cp - sp)*(cp + sp)*(-2 + 3*CT[4] + 3*ST[2] + 15*ST[4] - CT[2]*(1 + 30*ST[2])))/128.;
			break;
		case 29:
			return (cp*sqrtl(165/(2.*M_PI))*st*(2 - 7*ST[2] + 21*(5*CT[4] + ST[4] + CT[2]*(1 - 10*ST[2]))))/256.;
			break;
		case 30:
			return (ct*sqrtl(11/M_PI)*(30 + 35*CT[2] + 63*CT[4] - 105*(1 + 6*CT[2])*ST[2] + 315*ST[4]))/256.;
			break;
		case 31:
			return -(cp*sqrtl(165/(2.*M_PI))*st*(2 - 7*ST[2] + 21*(5*CT[4] + ST[4] + CT[2]*(1 - 10*ST[2]))))/256.;
			break;
		case 32:
			return -(ct*sqrtl(1155/(2.*M_PI))*(cp - sp)*(cp + sp)*(-2 + 3*CT[4] + 3*ST[2] + 15*ST[4] - CT[2]*(1 + 30*ST[2])))/128.;
			break;
		case 33:
			return (cp*sqrtl(385/M_PI)*(CP[2] - 3*SP[2])*st*(-6 + 45*CT[4] + 13*ST[2] + 9*ST[4] - 3*CT[2]*(13 + 30*ST[2])))/512.;
			break;
		case 34:
			return (3*ct*sqrtl(385/(2.*M_PI))*(CP[4] - 6*CP[2]*SP[2] + SP[4])*(2 + CT[4] + 9*ST[2] + 5*ST[4] - CT[2]*(3 + 10*ST[2])))/256.;
			break;
		case 35:
			return (-3*cp*sqrtl(77/M_PI)*(CP[4] - 10*CP[2]*SP[2] + 5*SP[4])*st*(10 + 5*CT[4] + 5*ST[2] + ST[4] - 5*CT[2]*(3 + 2*ST[2])))/512.;
			break;
	
	//l=6
		case 36:
			return -(sqrtl(3003/M_PI)*(CP[6] - 15*CP[4]*SP[2] + 15*CP[2]*SP[4] - SP[6])*(CT[6] - 3*CT[4]*(2 + 5*ST[2]) - (1 + ST[2])*(10 + 5*ST[2] + ST[4]) + 3*CT[2]*(5 + 12*ST[2] + 5*ST[4])))/2048.;
			break;
		case 37:
			return (3*cp*ct*sqrtl(1001/M_PI)*(CP[4] - 10*CP[2]*SP[2] + 5*SP[4])*st*(5 + 3*CT[4] + 8*ST[2] + 3*ST[4] - 2*CT[2]*(4 + 5*ST[2])))/512.;
			break;
		case 38:
			return (3*sqrtl(91/(2.*M_PI))*(CP[4] - 6*CP[2]*SP[2] + SP[4])*(10 + 11*CT[6] - 5*ST[2] - 26*ST[4] - 11*ST[6] - CT[4]*(26 + 165*ST[2]) + CT[2]*(5 + 156*ST[2] + 165*ST[4])))/1024.;
			break;
		case 39:
			return -(cp*ct*sqrtl(1365/M_PI)*(CP[2] - 3*SP[2])*st*(-9 + 33*CT[4] + 24*ST[2] + 33*ST[4] - 2*CT[2]*(12 + 55*ST[2])))/512.;
			break;
		case 40:
			return -(sqrtl(1365/M_PI)*(cp - sp)*(cp + sp)*(-10 + 33*CT[6] + 17*ST[2] - 6*ST[4] - 33*ST[6] - 3*CT[4]*(2 + 165*ST[2]) + CT[2]*(-17 + 36*ST[2] + 495*ST[4])))/2048.;
			break;
		case 41:
			return (cp*ct*sqrtl(273/(2.*M_PI))*st*(5 + 24*CT[2] + 99*CT[4] - 6*(4 + 55*CT[2])*ST[2] + 99*ST[4]))/256.;
			break;
		case 42:
			return (sqrtl(13/M_PI)*(50 - 105*ST[2] + 21*(11*CT[6] + 6*ST[4] - 11*ST[6] + CT[4]*(6 - 165*ST[2]) + CT[2]*(5 - 36*ST[2] + 165*ST[4]))))/1024.;
			break;
		case 43:
			return -(cp*ct*sqrtl(273/(2.*M_PI))*st*(5 + 24*CT[2] + 99*CT[4] - 6*(4 + 55*CT[2])*ST[2] + 99*ST[4]))/256.;
			break;
		case 44:
			return -(sqrtl(1365/M_PI)*(cp - sp)*(cp + sp)*(-10 + 33*CT[6] + 17*ST[2] - 6*ST[4] - 33*ST[6] - 3*CT[4]*(2 + 165*ST[2]) + CT[2]*(-17 + 36*ST[2] + 495*ST[4])))/2048.;
			break;
		case 45:
			return (cp*ct*sqrtl(1365/M_PI)*(CP[2] - 3*SP[2])*st*(-9 + 33*CT[4] + 24*ST[2] + 33*ST[4] - 2*CT[2]*(12 + 55*ST[2])))/512.;
			break;
		case 46:
			return (3*sqrtl(91/(2.*M_PI))*(CP[4] - 6*CP[2]*SP[2] + SP[4])*(10 + 11*CT[6] - 5*ST[2] - 26*ST[4] - 11*ST[6] - CT[4]*(26 + 165*ST[2]) + CT[2]*(5 + 156*ST[2] + 165*ST[4])))/1024.;
			break;
		case 47:
			return (-3*cp*ct*sqrtl(1001/M_PI)*(CP[4] - 10*CP[2]*SP[2] + 5*SP[4])*st*(5 + 3*CT[4] + 8*ST[2] + 3*ST[4] - 2*CT[2]*(4 + 5*ST[2])))/512.;
			break;
		case 48:
			return -(sqrtl(3003/M_PI)*(CP[6] - 15*CP[4]*SP[2] + 15*CP[2]*SP[4] - SP[6])*(CT[6] - 3*CT[4]*(2 + 5*ST[2]) - (1 + ST[2])*(10 + 5*ST[2] + ST[4]) + 3*CT[2]*(5 + 12*ST[2] + 5*ST[4])))/2048.;
			break;
		default:
			return 0;
			break;
	}
}

gNum cSphY(gIdx lmIdx, gNum cp, gNum sp, gNum ct, gNum st){
	gNum SP[7]={1.};
	gNum ST[7]={1.};
	gNum CP[7]={1.};
	gNum CT[7]={1.};
	for(int i=1;i<=6;i++){
		SP[i]=sp*SP[i-1];
		ST[i]=st*ST[i-1];
		CP[i]=cp*CP[i-1];
		CT[i]=ct*CT[i-1];
	}

	switch(lmIdx)
	{

		//l=next
		
		case 0:
			return 0;
			break;
		//l=next
		
		case 1:
			return -(sqrtl(3/(2.*M_PI))*sp*st)/2.;
			break;
		case 2:
			return 0;
			break;
		case 3:
			return -(sqrtl(3/(2.*M_PI))*sp*st)/2.;
			break;
		//l=next

		case 4:
			return (cp*sqrtl(15/(2.*M_PI))*sp*(-1 + CT[2] - ST[2]))/4.;
			break;
		case 5:
			return -(ct*sqrtl(15/(2.*M_PI))*sp*st)/2.;
			break;
		case 6:
			return 0;
			break;
		case 7:
			return -(ct*sqrtl(15/(2.*M_PI))*sp*st)/2.;
			break;
		case 8:
			return (cp*sqrtl(15/(2.*M_PI))*sp*(1 - CT[2] + ST[2]))/4.;
			break;
		//l=next

		case 9:
			return (sqrtl(35/M_PI)*sp*(-3*CP[2] + SP[2])*st*(3 - 3*CT[2] + ST[2]))/32.;
			break;
		case 10:
			return (cp*ct*sqrtl(105/(2.*M_PI))*sp*(-1 + CT[2] - 3*ST[2]))/8.;
			break;
		case 11:
			return (sqrtl(21/M_PI)*sp*st*(-1 - 15*CT[2] + 5*ST[2]))/32.;
			break;
		case 12:
			return 0;
			break;
		case 13:
			return (sqrtl(21/M_PI)*sp*st*(-1 - 15*CT[2] + 5*ST[2]))/32.;
			break;
		case 14:
			return -(cp*ct*sqrtl(105/(2.*M_PI))*sp*(-1 + CT[2] - 3*ST[2]))/8.;
			break;
		case 15:
			return (sqrtl(35/M_PI)*sp*(-3*CP[2] + SP[2])*st*(3 - 3*CT[2] + ST[2]))/32.;
			break;
		//l=next

		case 16:
			return (-3*cp*sqrtl(35/(2.*M_PI))*(cp - sp)*sp*(cp + sp)*(3 + CT[4] + 4*ST[2] + ST[4] - 2*CT[2]*(2 + 3*ST[2])))/32.;
			break;
		case 17:
			return (-3*ct*sqrtl(35/M_PI)*sp*(-3*CP[2] + SP[2])*st*(-1 + CT[2] - ST[2]))/16.;
			break;
		case 18:
			return (3*cp*sqrtl(5/(2.*M_PI))*sp*(-3 + 7*CT[4] + 4*ST[2] + 7*ST[4] - 2*CT[2]*(2 + 21*ST[2])))/32.;
			break;
		case 19:
			return (-3*ct*sqrtl(5/M_PI)*sp*st*(1 + 7*CT[2] - 7*ST[2]))/16.;
			break;
		case 20:
			return 0;
			break;
		case 21:
			return (-3*ct*sqrtl(5/M_PI)*sp*st*(1 + 7*CT[2] - 7*ST[2]))/16.;
			break;
		case 22:
			return (3*cp*sqrtl(5/(2.*M_PI))*sp*(3 - 7*CT[4] - 4*ST[2] - 7*ST[4] + CT[2]*(4 + 42*ST[2])))/32.;
			break;
		case 23:
			return (-3*ct*sqrtl(35/M_PI)*sp*(-3*CP[2] + SP[2])*st*(-1 + CT[2] - ST[2]))/16.;
			break;
		case 24:
			return (3*cp*sqrtl(35/(2.*M_PI))*(cp - sp)*sp*(cp + sp)*(3 + CT[4] + 4*ST[2] + ST[4] - 2*CT[2]*(2 + 3*ST[2])))/32.;
			break;
		//l=next

		case 25:
			return (-3*sqrtl(77/M_PI)*sp*(5*CP[4] - 10*CP[2]*SP[2] + SP[4])*st*(10 + 5*CT[4] + 5*ST[2] + ST[4] - 5*CT[2]*(3 + 2*ST[2])))/512.;
			break;
		case 26:
			return (-3*cp*ct*sqrtl(385/(2.*M_PI))*(cp - sp)*sp*(cp + sp)*(2 + CT[4] + 9*ST[2] + 5*ST[4] - CT[2]*(3 + 10*ST[2])))/64.;
			break;
		case 27:
			return -(sqrtl(385/M_PI)*sp*(-3*CP[2] + SP[2])*st*(-6 + 45*CT[4] + 13*ST[2] + 9*ST[4] - 3*CT[2]*(13 + 30*ST[2])))/512.;
			break;
		case 28:
			return (cp*ct*sqrtl(1155/(2.*M_PI))*sp*(-2 + 3*CT[4] + 3*ST[2] + 15*ST[4] - CT[2]*(1 + 30*ST[2])))/64.;
			break;
		case 29:
			return -(sqrtl(165/(2.*M_PI))*sp*st*(2 - 7*ST[2] + 21*(5*CT[4] + ST[4] + CT[2]*(1 - 10*ST[2]))))/256.;
			break;
		case 30:
			return 0;
			break;
		case 31:
			return -(sqrtl(165/(2.*M_PI))*sp*st*(2 - 7*ST[2] + 21*(5*CT[4] + ST[4] + CT[2]*(1 - 10*ST[2]))))/256.;
			break;
		case 32:
			return (cp*ct*sqrtl(1155/(2.*M_PI))*sp*(2 - 3*CT[4] + CT[2]*(1 + 30*ST[2]) - 3*(ST[2] + 5*ST[4])))/64.;
			break;
		case 33:
			return -(sqrtl(385/M_PI)*sp*(-3*CP[2] + SP[2])*st*(-6 + 45*CT[4] + 13*ST[2] + 9*ST[4] - 3*CT[2]*(13 + 30*ST[2])))/512.;
			break;
		case 34:
			return (3*cp*ct*sqrtl(385/(2.*M_PI))*(cp - sp)*sp*(cp + sp)*(2 + CT[4] + 9*ST[2] + 5*ST[4] - CT[2]*(3 + 10*ST[2])))/64.;
			break;
		case 35:
			return (-3*sqrtl(77/M_PI)*sp*(5*CP[4] - 10*CP[2]*SP[2] + SP[4])*st*(10 + 5*CT[4] + 5*ST[2] + ST[4] - 5*CT[2]*(3 + 2*ST[2])))/512.;
			break;
		//l=next

		case 36:
			return (cp*sqrtl(3003/M_PI)*sp*(3*CP[4] - 10*CP[2]*SP[2] + 3*SP[4])*(CT[6] - 3*CT[4]*(2 + 5*ST[2]) - (1 + ST[2])*(10 + 5*ST[2] + ST[4]) + 3*CT[2]*(5 + 12*ST[2] + 5*ST[4])))/1024.;
			break;
		case 37:
			return (-3*ct*sqrtl(1001/M_PI)*sp*(5*CP[4] - 10*CP[2]*SP[2] + SP[4])*st*(5 + 3*CT[4] + 8*ST[2] + 3*ST[4] - 2*CT[2]*(4 + 5*ST[2])))/512.;
			break;
		case 38:
			return (-3*cp*sqrtl(91/(2.*M_PI))*(cp - sp)*sp*(cp + sp)*(10 + 11*CT[6] - 5*ST[2] - 26*ST[4] - 11*ST[6] - CT[4]*(26 + 165*ST[2]) + CT[2]*(5 + 156*ST[2] + 165*ST[4])))/256.;
			break;
		case 39:
			return -(ct*sqrtl(1365/M_PI)*sp*(-3*CP[2] + SP[2])*st*(-9 + 33*CT[4] + 24*ST[2] + 33*ST[4] - 2*CT[2]*(12 + 55*ST[2])))/512.;
			break;
		case 40:
			return (cp*sqrtl(1365/M_PI)*sp*(-10 + 33*CT[6] + 17*ST[2] - 6*ST[4] - 33*ST[6] - 3*CT[4]*(2 + 165*ST[2]) + CT[2]*(-17 + 36*ST[2] + 495*ST[4])))/1024.;
			break;
		case 41:
			return -(ct*sqrtl(273/(2.*M_PI))*sp*st*(5 + 24*CT[2] + 99*CT[4] - 6*(4 + 55*CT[2])*ST[2] + 99*ST[4]))/256.;
			break;
		case 42:
			return 0;
			break;
		case 43:
			return -(ct*sqrtl(273/(2.*M_PI))*sp*st*(5 + 24*CT[2] + 99*CT[4] - 6*(4 + 55*CT[2])*ST[2] + 99*ST[4]))/256.;
			break;
		case 44:
			return (cp*sqrtl(1365/M_PI)*sp*(10 - 33*CT[6] - 17*ST[2] + 6*ST[4] + 33*ST[6] + CT[4]*(6 + 495*ST[2]) + CT[2]*(17 - 9*ST[2]*(4 + 55*ST[2]))))/1024.;
			break;
		case 45:
			return -(ct*sqrtl(1365/M_PI)*sp*(-3*CP[2] + SP[2])*st*(-9 + 33*CT[4] + 24*ST[2] + 33*ST[4] - 2*CT[2]*(12 + 55*ST[2])))/512.;
			break;
		case 46:
			return (3*cp*sqrtl(91/(2.*M_PI))*(cp - sp)*sp*(cp + sp)*(10 + 11*CT[6] - 5*ST[2] - 26*ST[4] - 11*ST[6] - CT[4]*(26 + 165*ST[2]) + CT[2]*(5 + 156*ST[2] + 165*ST[4])))/256.;
			break;
		case 47:
			return (-3*ct*sqrtl(1001/M_PI)*sp*(5*CP[4] - 10*CP[2]*SP[2] + SP[4])*st*(5 + 3*CT[4] + 8*ST[2] + 3*ST[4] - 2*CT[2]*(4 + 5*ST[2])))/512.;
			break;
		case 48:
			return -(cp*sqrtl(3003/M_PI)*sp*(3*CP[4] - 10*CP[2]*SP[2] + 3*SP[4])*(CT[6] - 3*CT[4]*(2 + 5*ST[2]) - (1 + ST[2])*(10 + 5*ST[2] + ST[4]) + 3*CT[2]*(5 + 12*ST[2] + 5*ST[4])))/1024.;
			break;
		default:
			return 0;
			break;
	}
}


//Old Spherical Harmonic functions

/*

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

*/