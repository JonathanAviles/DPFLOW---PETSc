#include "pf.h"

#undef __FUNCT__
#define __FUNCT__ "PrintDPF"
PetscErrorCode PrintDPF(PFDATA *pfdata){
	PetscErrorCode ierr;
	//PetscInt i, j;

	PetscFunctionBegin;

	ierr = PetscPrintf(PETSC_COMM_SELF,"**** Power flow dist case ****\n\n");CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_SELF,"Base power = %lf, nbus = %d, ngen = %d, nwye = %d, ndelta = %d, nbranch = %d\n",\
		pfdata->sbase,pfdata->nbus,pfdata->ngen,pfdata->nwye,pfdata->ndelta,pfdata->nbranch);CHKERRQ(ierr);
/*
	ierr = PetscPrintf(PETSC_COMM_SELF,"Bus data:\n");CHKERRQ(ierr);
	for (i=0; i<pfdata->nbus; i++){
		ierr = PetscPrintf(PETSC_COMM_SELF,"ID=%d\t",i);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_SELF,"base=%lf\t",pfdata->bus[i].basekV);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_SELF,"type=%d\t",pfdata->bus[i].ide);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_SELF,"gs=%lf\t",pfdata->bus[i].gl);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_SELF,"bs=%lf\t",pfdata->bus[i].bl);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_SELF,"vm=%lf\t",pfdata->bus[i].vm);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_SELF,"va=%lf\t",pfdata->bus[i].va);CHKERRQ(ierr);

		ierr = PetscPrintf(PETSC_COMM_SELF,"ngen=%d\t",pfdata->bus[i].ngen);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_SELF,"Gen IDs: ");CHKERRQ(ierr);
		if (pfdata->bus[i].ngen == 0) {
			ierr = PetscPrintf(PETSC_COMM_SELF,"\t");CHKERRQ(ierr);
		} else {
			for (j=0; j<pfdata->bus[i].ngen; j++){
				ierr = PetscPrintf(PETSC_COMM_SELF,"%d\t",pfdata->bus[i].gidx[j]);CHKERRQ(ierr);
			}
		}

		ierr = PetscPrintf(PETSC_COMM_SELF,"nload=%d\t",pfdata->bus[i].nload);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_SELF,"Wye-load IDs: ");CHKERRQ(ierr);
		if (pfdata->bus[i].nload == 0) {
			ierr = PetscPrintf(PETSC_COMM_SELF,"\t");CHKERRQ(ierr);
		} else {
			for (j=0; j<pfdata->bus[i].nload; j++){
				ierr = PetscPrintf(PETSC_COMM_SELF,"%d\t",pfdata->bus[i].lidx[j]);CHKERRQ(ierr);
			}
		}

		ierr = PetscPrintf(PETSC_COMM_SELF,"ndload=%d\t",pfdata->bus[i].ndload);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_SELF,"Delta-load IDs: ");CHKERRQ(ierr);
		if (pfdata->bus[i].ndload == 0) {
			ierr = PetscPrintf(PETSC_COMM_SELF,"\t");CHKERRQ(ierr);
		} else {
			for (j=0; j<pfdata->bus[i].ndload; j++){
				ierr = PetscPrintf(PETSC_COMM_SELF,"%d\t",pfdata->bus[i].dlidx[j]);CHKERRQ(ierr);
			}
		}

		ierr = PetscPrintf(PETSC_COMM_SELF,"\n");CHKERRQ(ierr);
	}
	ierr = PetscPrintf(PETSC_COMM_SELF,"\n");CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_SELF,"Gen data:\n");CHKERRQ(ierr);
	for (i=0; i<pfdata->ngen; i++){
		ierr = PetscPrintf(PETSC_COMM_SELF,"ID=%d\t",i);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_SELF,"Pg=%lf\t",pfdata->gen[i].pg);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_SELF,"Qg=%lf\n",pfdata->gen[i].qg);CHKERRQ(ierr);
	}
	ierr = PetscPrintf(PETSC_COMM_SELF,"\n");CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_SELF,"Wye-load data:\n");CHKERRQ(ierr);
	for (i=0; i<pfdata->nwye; i++){
		ierr = PetscPrintf(PETSC_COMM_SELF,"ID=%d\t",i);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_SELF,"pl=%lf\t",pfdata->wye[i].pl);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_SELF,"ql=%lf\t",pfdata->wye[i].ql);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_SELF,"ip=%lf\t",pfdata->wye[i].ip);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_SELF,"iq=%lf\n",pfdata->wye[i].iq);CHKERRQ(ierr);
	}
	ierr = PetscPrintf(PETSC_COMM_SELF,"\n");CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_SELF,"Delta-load data:\n");CHKERRQ(ierr);
	for (i=0; i<pfdata->ndelta; i++){
		ierr = PetscPrintf(PETSC_COMM_SELF,"ID=%d\t",i);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_SELF,"pl=%lf\t",pfdata->delta[i].pl);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_SELF,"ql=%lf\t",pfdata->delta[i].ql);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_SELF,"ip=%lf\t",pfdata->delta[i].ip);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_SELF,"iq=%lf\n",pfdata->delta[i].iq);CHKERRQ(ierr);
	}
	ierr = PetscPrintf(PETSC_COMM_SELF,"\n");CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_SELF,"Branch data:\n");CHKERRQ(ierr);
	for (i=0; i<pfdata->nbranch; i++){
		if (pfdata->branch[i].status){
			ierr = PetscPrintf(PETSC_COMM_SELF,"#=%d\t",i);CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_SELF,"internal_i=%d\t",pfdata->branch[i].internal_i);CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_SELF,"internal_j=%d\t",pfdata->branch[i].internal_j);CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_SELF,"yff=%lf+i*%lf\t",pfdata->branch[i].yff[0],pfdata->branch[i].yff[1]);CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_SELF,"yft=%lf+i*%lf\t",pfdata->branch[i].yft[0],pfdata->branch[i].yft[1]);CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_SELF,"ytf=%lf+i*%lf\t",pfdata->branch[i].ytf[0],pfdata->branch[i].ytf[1]);CHKERRQ(ierr);
			ierr = PetscPrintf(PETSC_COMM_SELF,"ytt=%lf+i*%lf\n",pfdata->branch[i].ytt[0],pfdata->branch[i].ytt[1]);CHKERRQ(ierr);
		}
	}
	ierr = PetscPrintf(PETSC_COMM_SELF,"\n");CHKERRQ(ierr);
*/
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DPFReadMatPowerData"
PetscErrorCode DPFReadMatPowerData(PFDATA *pf, char *fileName)
{
	FILE							*fp;
	char							fileLine[MAXCHARSPERLINE];
 	PetscErrorCode		ierr;
 	VERTEXDATA				equivBus;
	WYELOAD						equivWyeLoad, equivWyeLoad_ord;
	DELTALOAD					equivDeltaLoad;
	GEN								equivGen, equivGen_ord;
	EDGEDATA					equivBranch;
	PetscInt					extBusID, bus_Type, gen_busID, gen_status;
	PetscScalar				bus_Vmag[3], bus_Vang[3], baseVoltage, gen_P[3], gen_Q[3];
	PetscInt					wye_busID, wye_status, delta_busID, delta_status;
	PetscScalar				wye_ctPow[6], wye_ctCur[6],wye_ctImp[6], delta_ctPow[6], delta_ctCur[6], delta_ctImp[6];
  PetscInt					line_fbus, line_tbus, line_status, trafo_fbus, trafo_tbus, trafoConn, trafo_status;
  PetscScalar				Rbr[6], Xbr[6], Bbr[6], Rtr, Xtr, primTap, secTap;
	PetscInt					numBus, numGen, numWyeLoad, numDeltaLoad, numLineBranch, numTrafoBranch;
	PetscInt					i=0, j=0, j1=0, j2=0, j3=0, j4=0, j5=0, k=0, m=0, numNonConPhases;
	PetscInt					lineCounter=0, maxBusID=-1, offset;
	PetscInt					bus_start_line=-1,		bus_end_line=-1;
	PetscInt					gen_start_line=-1,		gen_end_line=-1;
	PetscInt					wye_start_line=-1,		wye_end_line=-1;
	PetscInt					delta_start_line=-1,	delta_end_line=-1;
	PetscInt					line_start_line=-1,		line_end_line=-1;
	PetscInt					trafo_start_line=-1,	trafo_end_line=-1;
	PetscInt					*ext2int_busMap;
	PetscInt					noPhaseIndexes[3], indexFromTo[6];
	PetscScalar				Ybr_real[9], Ybr_imag[9], Zbr_real[9], Zbr_imag[9], Bbr_3x3[9];
	PetscScalar				auxZ2x2_real[4], auxZ2x2_imag[4], auxY2x2_real[4], auxY2x2_imag[4];
	PetscScalar				basePower, ytr_pu_real, ytr_pu_imag;
	PetscScalar				Ypp_real[9], Ypp_imag[9], Yps_real[9], Yps_imag[9], Ysp_real[9], Ysp_imag[9], Yss_real[9], Yss_imag[9];
	PetscScalar				lineCond_val[15], lineSusc_val[15];
	const PetscScalar	YI[9]    = {1,0,0,0,1,0,0,0,1};
	const PetscScalar	YII[9]   = {2/3.0,-1/3.0,-1/3.0,-1/3.0,2/3.0,-1/3.0,-1/3.0,-1/3.0,2/3.0};
	const PetscScalar	YIII[9]  = {-1/sqrt(3),1/sqrt(3),0,0,-1/sqrt(3),1/sqrt(3),1/sqrt(3),0,-1/sqrt(3)};
	const PetscScalar	YIIIt[9] = {-1/sqrt(3),0,1/sqrt(3),1/sqrt(3),-1/sqrt(3),0,0,1/sqrt(3),-1/sqrt(3)};

	/* Initialize PETSc function */
	PetscFunctionBegin;

	/****************************************************** Initial calculations ************************************************************
	 - Base power is read from data file.
	 - Number of buses, generators, wye-loads, delta-loads, line branches and transformer branches are determined.
	*****************************************************************************************************************************************/

	/* Read MATPOWER casefile */
	fp = fopen(fileName,"r");
	if (!fp)	SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Can't open data file %s",fileName); /*Check if it is a valid file*/

	/* File is read line by line to determine where the matrices caseObject.bus, caseObject.gen, caseObject.yload,
		 caseObject.dload, caseObject.lineBranch and caseObject.trafoBranch start and end */
	while (fgets(fileLine,MAXCHARSPERLINE,fp)) {

		/* Setting pf->sbase according to the case file */
		if (strstr(fileLine,"caseObject.baseMVA")){
			sscanf(fileLine,"%*[^0123456789]%lf",&basePower);
			pf->sbase = basePower;
		}

		/* Line where matrix starts is determined */
		if (strstr(fileLine,"caseObject.bus") && bus_start_line == -1)						bus_start_line = lineCounter+1; /* Bus data starts from next line */
		if (strstr(fileLine,"caseObject.gen") && gen_start_line == -1)						gen_start_line = lineCounter+1; /* Bus data starts from next line */
		if (strstr(fileLine,"caseObject.yload") && wye_start_line == -1)					wye_start_line = lineCounter+1; /* Bus data starts from next line */
		if (strstr(fileLine,"caseObject.dload") && delta_start_line == -1)				delta_start_line = lineCounter+1; /* Bus data starts from next line */
		if (strstr(fileLine,"caseObject.lineBranch") && line_start_line == -1)		line_start_line = lineCounter+1; /* Bus data starts from next line */
		if (strstr(fileLine,"caseObject.trafoBranch") && trafo_start_line == -1)	trafo_start_line = lineCounter+1; /* Bus data starts from next line */

		/* Line where matrix finishes is determined */
		if (strstr(fileLine,"];")) {
			if (bus_start_line != -1 && bus_end_line == -1)			bus_end_line = lineCounter;
			if (gen_start_line != -1 && gen_end_line == -1)			gen_end_line = lineCounter;
			if (wye_start_line != -1 && wye_end_line == -1)			wye_end_line = lineCounter;
			if (delta_start_line != -1 && delta_end_line == -1)	delta_end_line = lineCounter;
			if (line_start_line != -1 && line_end_line == -1)		line_end_line = lineCounter;
			if (trafo_start_line != -1 && trafo_end_line == -1)	trafo_end_line = lineCounter;
		}

		/* Determine the total number of buses */
		if (bus_start_line != -1 && lineCounter >= bus_start_line && bus_end_line == -1) {
			sscanf(fileLine,"%d",&extBusID);
			if (extBusID > maxBusID) maxBusID = extBusID;
		}
		lineCounter++;
	} /* End of while (file reading) */
	fclose(fp);

	/* Numbers of buses, generators, yloads, dloads, line branches and transformer branches
	   are determined from the line numbers previously obtained */
	numBus = bus_end_line - bus_start_line;
	numGen = gen_end_line - gen_start_line;
	numWyeLoad = wye_end_line - wye_start_line;
	numDeltaLoad = delta_end_line - delta_start_line;
	numLineBranch = line_end_line - line_start_line;
	numTrafoBranch = trafo_end_line - trafo_start_line;

		/* Number of elements is printed */
	ierr = PetscPrintf(PETSC_COMM_SELF,\
		"Base power = %lf, numbus = %d, numgen = %d, numyl = %d, numdl = %d, numlbr = %d, numtbr = %d\n\n",\
		basePower, numBus, numGen, numWyeLoad, numDeltaLoad, numLineBranch, numTrafoBranch);CHKERRQ(ierr);

	/* Number of buses, generators and loads are assigned to the power flow structure */
	pf->nbus    = 3*numBus;
	pf->ngen    = 3*numGen;
	pf->nwye		= 3*numWyeLoad;
	pf->ndelta	= 3*numDeltaLoad;
	pf->nbranch = 15*(numLineBranch + numTrafoBranch) + 3*numDeltaLoad;

	/* Using the total number of elements, memory is allocated */
	ierr = PetscCalloc1(pf->nbus,&pf->bus);CHKERRQ(ierr);
	ierr = PetscCalloc1(pf->ngen,&pf->gen);CHKERRQ(ierr);
	ierr = PetscCalloc1(pf->nwye,&pf->wye);CHKERRQ(ierr);
	ierr = PetscCalloc1(pf->ndelta,&pf->delta);CHKERRQ(ierr);
	ierr = PetscCalloc1(pf->nbranch,&pf->branch);CHKERRQ(ierr);

	/* Internal pointers are defined for the function (bus, generator and load) */
	equivBus = pf->bus; equivGen = pf->gen; equivWyeLoad = pf->wye; equivDeltaLoad = pf->delta; equivBranch = pf->branch;

	/* Memory allocation for index mapping vector (external to 3-phase internal)*/
	ierr = PetscMalloc1(maxBusID+1,&ext2int_busMap);CHKERRQ(ierr);
	for (i=0; i <= maxBusID; i++) ext2int_busMap[i] = -1;

	/* Initialize number of generators and loads connected to each bus */
	for(i=0; i < pf->nbus; i++) {
		equivBus[i].ngen = equivBus[i].nload = equivBus[i].ndload = equivBus[i].isConnBus =	0;
		equivBus[i].gl = equivBus[i].bl = 0.0;
	}

	offset = pf->nbranch-3*numDeltaLoad;

	/*********************************************************** Data lecture ***************************************************************
	 - Memory is allocated for data from dist case file
	 - Data is read from dist case file
	*****************************************************************************************************************************************/

	/* Reading data from case file */
	fp = fopen(fileName,"r");
	for (i=0;i<lineCounter;i++){
		/* File is read line by line again */
		fgets(fileLine,MAXCHARSPERLINE,fp);

		/************* Bus data *************/
		if ((i >= bus_start_line) && (i < bus_end_line)){
			sscanf(fileLine,"%d %d %lf %lf %lf %lf %lf %lf %lf",&extBusID,&bus_Type,\
				&bus_Vmag[0],&bus_Vmag[1],&bus_Vmag[2],&bus_Vang[0],&bus_Vang[1],&bus_Vang[2],&baseVoltage);
			j1++;
			ext2int_busMap[extBusID] = j1;
			for (j=0;j<3;j++){
				equivBus[3*j1-3+j].ide = bus_Type;
				equivBus[3*j1-3+j].vm = bus_Vmag[j];
				equivBus[3*j1-3+j].va = bus_Vang[j];
				equivBus[3*j1-3+j].basekV = baseVoltage;
			}
		}

		/********** Generator data **********/
		if ((i >= gen_start_line) && (i < gen_end_line)){
			sscanf(fileLine,"%d %lf %lf %lf %lf %lf %lf %d",&gen_busID,\
				&gen_P[0],&gen_P[1],&gen_P[2],&gen_Q[0],&gen_Q[1],&gen_Q[2],&gen_status);
			j2++;
			if (gen_status!=0){
				for (j=0;j<3;j++){
					equivGen[3*j2-3+j].pg = gen_P[j];
					equivGen[3*j2-3+j].qg = gen_Q[j];
					equivBus[3*ext2int_busMap[gen_busID]-3+j].gidx[0]=3*j2-3+j;
					equivBus[3*ext2int_busMap[gen_busID]-3+j].ngen=1;
				}
			}
  	}

		/********** Wye-load data ***********/
		if ((i >= wye_start_line) && (i < wye_end_line)){
			sscanf(fileLine,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d",&wye_busID,\
				&wye_ctPow[0],&wye_ctPow[1],&wye_ctPow[2],&wye_ctPow[3],&wye_ctPow[4],&wye_ctPow[5],\
				&wye_ctCur[0],&wye_ctCur[1],&wye_ctCur[2],&wye_ctCur[3],&wye_ctCur[4],&wye_ctCur[5],\
				&wye_ctImp[0],&wye_ctImp[1],&wye_ctImp[2],&wye_ctImp[3],&wye_ctImp[4],&wye_ctImp[5],&wye_status);
			j3++;
			if(wye_status!=0){
				for(j=0; j<3; j++){
					equivWyeLoad[3*j3-3+j].pl = wye_ctPow[j];
					equivWyeLoad[3*j3-3+j].ql = wye_ctPow[j+3];
					equivWyeLoad[3*j3-3+j].ip = wye_ctCur[j];
					equivWyeLoad[3*j3-3+j].iq = wye_ctCur[j+3];
					equivBus[3*ext2int_busMap[wye_busID]-3+j].lidx[0] = 3*j3-3+j;
					equivBus[3*ext2int_busMap[wye_busID]-3+j].nload   = 1;
					equivBus[3*ext2int_busMap[wye_busID]-3+j].gl     += wye_ctImp[j];
					equivBus[3*ext2int_busMap[wye_busID]-3+j].bl     -= wye_ctImp[j+3];
				}
			}
		} /* End of wye-load data*/

		/********* Delta-load data **********/
		if ((i >= delta_start_line) && (i < delta_end_line)){
			sscanf(fileLine,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d",&delta_busID,\
				&delta_ctPow[0],&delta_ctPow[1],&delta_ctPow[2],&delta_ctPow[3],&delta_ctPow[4],&delta_ctPow[5],\
				&delta_ctCur[0],&delta_ctCur[1],&delta_ctCur[2],&delta_ctCur[3],&delta_ctCur[4],&delta_ctCur[5],\
				&delta_ctImp[0],&delta_ctImp[1],&delta_ctImp[2],&delta_ctImp[3],&delta_ctImp[4],&delta_ctImp[5],&delta_status);
			j4++;
			if(delta_status!=0){
				for(j=0; j<3; j++){
					equivDeltaLoad[3*j4-3+j].pl = delta_ctPow[j];
					equivDeltaLoad[3*j4-3+j].ql = delta_ctPow[j+3];
					equivDeltaLoad[3*j4-3+j].ip = delta_ctCur[j];
					equivDeltaLoad[3*j4-3+j].iq = delta_ctCur[j+3];
					equivBus[3*ext2int_busMap[delta_busID]-3+j].dlidx[0] = 3*j4-3+j;
					equivBus[3*ext2int_busMap[delta_busID]-3+j].ndload   = 1;

					equivBranch[offset+3*j4-3+j].status = 1;
					equivBranch[offset+3*j4-3+j].yff[0] = equivBranch[offset+3*j4-3+j].ytt[0] = delta_ctImp[j  ]/(3*basePower);
					equivBranch[offset+3*j4-3+j].yff[1] = equivBranch[offset+3*j4-3+j].ytt[1] = -delta_ctImp[j+3]/(3*basePower);
					equivBranch[offset+3*j4-3+j].yft[0] = equivBranch[offset+3*j4-3+j].ytf[0] = -delta_ctImp[j  ]/(3*basePower);
					equivBranch[offset+3*j4-3+j].yft[1] = equivBranch[offset+3*j4-3+j].ytf[1] = delta_ctImp[j+3]/(3*basePower);
				}

				equivBranch[offset+3*j4-3  ].internal_i = 3*ext2int_busMap[delta_busID]-3;
				equivBranch[offset+3*j4-3+1].internal_i = 3*ext2int_busMap[delta_busID]-2;
				equivBranch[offset+3*j4-3+2].internal_i = 3*ext2int_busMap[delta_busID]-1;

				equivBranch[offset+3*j4-3  ].internal_j = 3*ext2int_busMap[delta_busID]-2;
				equivBranch[offset+3*j4-3+2].internal_j = 3*ext2int_busMap[delta_busID]-3;
				equivBranch[offset+3*j4-3+1].internal_j = 3*ext2int_busMap[delta_busID]-1;
			} /* End of delta_status conditional */
		} /* End of delta-load data*/

		/********* Line-branch data *********/
		if ((i >= line_start_line) && (i < line_end_line)){
			sscanf(fileLine,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d",\
				&line_fbus,&line_tbus,&Rbr[0],&Rbr[1],&Rbr[2],&Rbr[3],&Rbr[4],&Rbr[5],&Xbr[0],&Xbr[1],&Xbr[2],&Xbr[3],\
				&Xbr[4],&Xbr[5],&Bbr[0],&Bbr[1],&Bbr[2],&Bbr[3],&Bbr[4],&Bbr[5],&line_status);

			if(line_status!=0){
				/* Computation of the 3x3 impedance matrix of each branch */
				Zbr_real[0] = Rbr[0]; Zbr_real[1] = Rbr[1]; Zbr_real[2] = Rbr[2];
				Zbr_real[3] = Rbr[1]; Zbr_real[4] = Rbr[3]; Zbr_real[5] = Rbr[4];
				Zbr_real[6] = Rbr[2]; Zbr_real[7] = Rbr[4]; Zbr_real[8] = Rbr[5];

				Zbr_imag[0] = Xbr[0]; Zbr_imag[1] = Xbr[1]; Zbr_imag[2] = Xbr[2];
				Zbr_imag[3] = Xbr[1]; Zbr_imag[4] = Xbr[3]; Zbr_imag[5] = Xbr[4];
				Zbr_imag[6] = Xbr[2]; Zbr_imag[7] = Xbr[4]; Zbr_imag[8] = Xbr[5];

				/* Computation of the 3x3 susceptance matrix of each branch */
				Bbr_3x3[0]  = Bbr[0]; Bbr_3x3[1]  = Bbr[1]; Bbr_3x3[2]  = Bbr[2];
				Bbr_3x3[3]  = Bbr[1]; Bbr_3x3[4]  = Bbr[3]; Bbr_3x3[5]  = Bbr[4];
				Bbr_3x3[6]  = Bbr[2]; Bbr_3x3[7]  = Bbr[4]; Bbr_3x3[8]  = Bbr[5];

				/* Initialize admittance matrix of branch to zero */
				for(j=0; j<9; j++){
					 Ybr_real[j] = Ybr_imag[j] = 0;
				}

				/* Computation of admittance matrix of each branch */
				zeroRowColIndexes_SymMat3x3(&numNonConPhases, noPhaseIndexes, Zbr_real, Zbr_imag);

				switch (numNonConPhases) {
					case 0:
						/* All phases are connected */
						ierr = complexMatrixInverse_3x3(Ybr_real, Ybr_imag, Zbr_real, Zbr_imag);CHKERRQ(ierr);
						break;

					case 1:
						/* Phase A and B connected */
						if ((noPhaseIndexes[0]==0) && (noPhaseIndexes[1]==0)) {
							auxZ2x2_real[0] = Zbr_real[0]; auxZ2x2_real[1] = Zbr_real[1];
							auxZ2x2_real[2] = Zbr_real[3]; auxZ2x2_real[3] = Zbr_real[4];

							auxZ2x2_imag[0] = Zbr_imag[0]; auxZ2x2_imag[1] = Zbr_imag[1];
							auxZ2x2_imag[2] = Zbr_imag[3]; auxZ2x2_imag[3] = Zbr_imag[4];

							ierr = complexMatrixInverse_2x2(auxY2x2_real, auxY2x2_imag, auxZ2x2_real, auxZ2x2_imag);CHKERRQ(ierr);

							Ybr_real[0] = auxY2x2_real[0]; Ybr_real[1] = auxY2x2_real[1];
							Ybr_real[3] = auxY2x2_real[2]; Ybr_real[4] = auxY2x2_real[3];

							Ybr_imag[0] = auxY2x2_imag[0]; Ybr_imag[1] = auxY2x2_imag[1];
							Ybr_imag[3] = auxY2x2_imag[2]; Ybr_imag[4] = auxY2x2_imag[3];
						}
						/* Phase A and C connected */
						if ((noPhaseIndexes[0]==0) && (noPhaseIndexes[2]==0)) {
							auxZ2x2_real[0] = Zbr_real[0]; auxZ2x2_real[1] = Zbr_real[2];
							auxZ2x2_real[2] = Zbr_real[6]; auxZ2x2_real[3] = Zbr_real[8];

							auxZ2x2_imag[0] = Zbr_imag[0]; auxZ2x2_imag[1] = Zbr_imag[2];
							auxZ2x2_imag[2] = Zbr_imag[6]; auxZ2x2_imag[3] = Zbr_imag[8];

							ierr = complexMatrixInverse_2x2(auxY2x2_real, auxY2x2_imag, auxZ2x2_real, auxZ2x2_imag);CHKERRQ(ierr);

							Ybr_real[0] = auxY2x2_real[0]; Ybr_real[2] = auxY2x2_real[1];
							Ybr_real[6] = auxY2x2_real[2]; Ybr_real[8] = auxY2x2_real[3];

							Ybr_imag[0] = auxY2x2_imag[0]; Ybr_imag[2] = auxY2x2_imag[1];
							Ybr_imag[6] = auxY2x2_imag[2]; Ybr_imag[8] = auxY2x2_imag[3];
						}
						/* Phase B and C connected */
						if ((noPhaseIndexes[1]==0) && (noPhaseIndexes[2]==0)) {
							auxZ2x2_real[0] = Zbr_real[4]; auxZ2x2_real[1] = Zbr_real[5];
							auxZ2x2_real[2] = Zbr_real[7]; auxZ2x2_real[3] = Zbr_real[8];

							auxZ2x2_imag[0] = Zbr_imag[4]; auxZ2x2_imag[1] = Zbr_imag[5];
							auxZ2x2_imag[2] = Zbr_imag[7]; auxZ2x2_imag[3] = Zbr_imag[8];

							ierr = complexMatrixInverse_2x2(auxY2x2_real, auxY2x2_imag, auxZ2x2_real, auxZ2x2_imag);CHKERRQ(ierr);

							Ybr_real[4] = auxY2x2_real[0]; Ybr_real[5] = auxY2x2_real[1];
							Ybr_real[7] = auxY2x2_real[2]; Ybr_real[8] = auxY2x2_real[3];

							Ybr_imag[4] = auxY2x2_imag[0]; Ybr_imag[5] = auxY2x2_imag[1];
							Ybr_imag[7] = auxY2x2_imag[2]; Ybr_imag[8] = auxY2x2_imag[3];
						}
						break;

					case 2:
						/* Only phase A is connected */
						if (noPhaseIndexes[0]==0) {
							ierr = complexScalarInverse(Ybr_real, Ybr_imag, Zbr_real[0], Zbr_imag[0]);CHKERRQ(ierr);
						}
						/* Only phase B is connected */
						if (noPhaseIndexes[1]==0) {
							ierr = complexScalarInverse(Ybr_real+4, Ybr_imag+4, Zbr_real[4], Zbr_imag[4]);CHKERRQ(ierr);
						}
						/* Only phase C is connected */
						if (noPhaseIndexes[2]==0) {
							ierr = complexScalarInverse(Ybr_real+8, Ybr_imag+8, Zbr_real[8], Zbr_imag[8]);CHKERRQ(ierr);
						}
						break;
				}

				for (j=0; j<3; j++){
					indexFromTo[j]   = 3*ext2int_busMap[line_fbus]-3+j;
					indexFromTo[j+3] = 3*ext2int_busMap[line_tbus]-3+j;
					if (noPhaseIndexes[j] == 0){
						equivBus[indexFromTo[j]].isConnBus   = 1;
						equivBus[indexFromTo[j+3]].isConnBus = 1;
					}
					equivBus[indexFromTo[j]].bl   += basePower*(Bbr_3x3[3*j]+Bbr_3x3[3*j+1]+Bbr_3x3[3*j+2])/2.0;
					equivBus[indexFromTo[j+3]].bl += basePower*(Bbr_3x3[3*j]+Bbr_3x3[3*j+1]+Bbr_3x3[3*j+2])/2.0;
				}

				m = 0;
				for (j=0; j<5; j++){
					for (k=j+1; k<6; k++){
						equivBranch[15*j5+m].internal_i = indexFromTo[j];
						equivBranch[15*j5+m].internal_j = indexFromTo[k];
						m++;
					}
				}

				lineCond_val[0]  = - Ybr_real[1];		lineSusc_val[0]  = - Ybr_imag[1] - Bbr_3x3[1]/2;
				lineCond_val[1]  = - Ybr_real[2];		lineSusc_val[1]  = - Ybr_imag[2] - Bbr_3x3[2]/2;
				lineCond_val[2]  = + Ybr_real[0];		lineSusc_val[2]  = + Ybr_imag[0];
				lineCond_val[3]  = + Ybr_real[1];		lineSusc_val[3]  = + Ybr_imag[1];
				lineCond_val[4]  = + Ybr_real[2];		lineSusc_val[4]  = + Ybr_imag[2];
				lineCond_val[5]  = - Ybr_real[5];		lineSusc_val[5]  = - Ybr_imag[5] - Bbr_3x3[5]/2;
				lineCond_val[6]  = + Ybr_real[3];		lineSusc_val[6]  = + Ybr_imag[3];
				lineCond_val[7]  = + Ybr_real[4];		lineSusc_val[7]  = + Ybr_imag[4];
				lineCond_val[8]  = + Ybr_real[5];		lineSusc_val[8]  = + Ybr_imag[5];
				lineCond_val[9]  = + Ybr_real[6];		lineSusc_val[9]  = + Ybr_imag[6];
				lineCond_val[10] = + Ybr_real[7];		lineSusc_val[10] = + Ybr_imag[7];
				lineCond_val[11] = + Ybr_real[8];		lineSusc_val[11] = + Ybr_imag[8];
				lineCond_val[12] = - Ybr_real[1];		lineSusc_val[12] = - Ybr_imag[1] - Bbr_3x3[1]/2;
				lineCond_val[13] = - Ybr_real[2];		lineSusc_val[13] = - Ybr_imag[2] - Bbr_3x3[2]/2;
				lineCond_val[14] = - Ybr_real[5];		lineSusc_val[14] = - Ybr_imag[5] - Bbr_3x3[5]/2;

				for (j=0; j<15; j++){
					equivBranch[15*j5+j].status = 1;
					equivBranch[15*j5+j].yff[0] = equivBranch[15*j5+j].ytt[0] = lineCond_val[j];
					equivBranch[15*j5+j].yff[1] = equivBranch[15*j5+j].ytt[1]	= lineSusc_val[j];
					equivBranch[15*j5+j].yft[0] = equivBranch[15*j5+j].ytf[0] = -lineCond_val[j];
					equivBranch[15*j5+j].yft[1] = equivBranch[15*j5+j].ytf[1] = -lineSusc_val[j];
				}

			} /* End of line conditional */

			j5++;
		} /* End of line-branch data lecture */

		/***** Transformer-branch data ******/
		if ((i >= trafo_start_line) && (i < trafo_end_line)){
      sscanf(fileLine,"%d %d %lf %lf %d %lf %lf %d",&trafo_fbus,&trafo_tbus,&Rtr,&Xtr,&trafoConn,&primTap,&secTap,&trafo_status);

			if (trafo_status != 0){
				/* Computer short-circuit admittance of transformer */
				ierr = complexScalarInverse(&ytr_pu_real, &ytr_pu_imag, Rtr, Xtr);CHKERRQ(ierr);

				switch (trafoConn) {
					case 1: /* YNyn0 */
						for(j=0; j<9; j++){
							Ypp_real[j] = ytr_pu_real*YI[j]/(primTap*primTap);  Ypp_imag[j] = ytr_pu_imag*YI[j]/(primTap*primTap);  /* YPP = (ytr_pu)*YI  */
							Yps_real[j] = -ytr_pu_real*YI[j]/(primTap*secTap);  Yps_imag[j] = -ytr_pu_imag*YI[j]/(primTap*secTap);  /* YPS = -(ytr_pu)*YI */
							Ysp_real[j] = -ytr_pu_real*YI[j]/(primTap*secTap);  Ysp_imag[j] = -ytr_pu_imag*YI[j]/(primTap*secTap);  /* YSP = -(ytr_pu)*YI */
							Yss_real[j] = ytr_pu_real*YI[j]/(secTap*secTap);    Yss_imag[j] = ytr_pu_imag*YI[j]/(secTap*secTap);    /* YSS = (ytr_pu)*YI  */
						}
						break;
					case 4: /* Dd0 */
						for(j=0; j<9; j++){
							Ypp_real[j] = ytr_pu_real*YII[j]/(primTap*primTap); Ypp_imag[j] = ytr_pu_imag*YII[j]/(primTap*primTap); /* YPP = (ytr_pu)*YII  */
							Yps_real[j] = -ytr_pu_real*YII[j]/(primTap*secTap); Yps_imag[j] = -ytr_pu_imag*YII[j]/(primTap*secTap); /* YPS = -(ytr_pu)*YII */
							Ysp_real[j] = -ytr_pu_real*YII[j]/(primTap*secTap); Ysp_imag[j] = -ytr_pu_imag*YII[j]/(primTap*secTap); /* YSP = -(ytr_pu)*YII */
							Yss_real[j] = ytr_pu_real*YII[j]/(secTap*secTap);   Yss_imag[j] = ytr_pu_imag*YII[j]/(secTap*secTap);   /* YSS = (ytr_pu)*YII  */
						}
						break;
					case 5: /* YNd1 */
						for(j=0; j<9; j++){
							Ypp_real[j] = ytr_pu_real*YI[j]/(primTap*primTap);  Ypp_imag[j] = ytr_pu_imag*YI[j]/(primTap*primTap); /* YPP = (ytr_pu)*YI   */
							Yps_real[j] = ytr_pu_real*YIII[j]/(primTap*secTap); Yps_imag[j] = ytr_pu_imag*YIII[j]/(primTap*secTap); /* YPS = (ytr_pu)*YIII  */
							Ysp_real[j] = ytr_pu_real*YIIIt[j]/(primTap*secTap);Ysp_imag[j] = ytr_pu_imag*YIIIt[j]/(primTap*secTap);/* YSP = (ytr_pu)*YIIIt */
							Yss_real[j] = ytr_pu_real*YII[j]/(secTap*secTap);   Yss_imag[j] = ytr_pu_imag*YII[j]/(secTap*secTap);    /* YSS = (ytr_pu)*YII    */
						}
						break;
					case 7: /* Dyn1 */
						for(j=0; j<9; j++){
							Ypp_real[j] = ytr_pu_real*YII[j]/(primTap*primTap); Ypp_imag[j] = ytr_pu_imag*YII[j]/(primTap*primTap); /* YPP = (ytr_pu)*YII   */
							Yps_real[j] = ytr_pu_real*YIII[j]/(primTap*secTap); Yps_imag[j] = ytr_pu_imag*YIII[j]/(primTap*secTap); /* YPS = (ytr_pu)*YIII  */
							Ysp_real[j] = ytr_pu_real*YIIIt[j]/(primTap*secTap);Ysp_imag[j] = ytr_pu_imag*YIIIt[j]/(primTap*secTap);/* YSP = (ytr_pu)*YIIIt */
							Yss_real[j] = ytr_pu_real*YI[j]/(secTap*secTap);    Yss_imag[j] = ytr_pu_imag*YI[j]/(secTap*secTap);    /* YSS = (ytr_pu)*YI    */
						}
						break;
					default: /* YNyn0 */
						for(j=0; j<9; j++){
							Ypp_real[j] = ytr_pu_real*YI[j]/(primTap*primTap);  Ypp_imag[j] = ytr_pu_imag*YI[j]/(primTap*primTap);  /* YPP = (ytr_pu)*YI  */
							Yps_real[j] = -ytr_pu_real*YI[j]/(primTap*secTap);  Yps_imag[j] = -ytr_pu_imag*YI[j]/(primTap*secTap);  /* YPS = -(ytr_pu)*YI */
							Ysp_real[j] = -ytr_pu_real*YI[j]/(primTap*secTap);  Ysp_imag[j] = -ytr_pu_imag*YI[j]/(primTap*secTap);  /* YSP = -(ytr_pu)*YI */
							Yss_real[j] = ytr_pu_real*YI[j]/(secTap*secTap);    Yss_imag[j] = ytr_pu_imag*YI[j]/(secTap*secTap);    /* YSS = (ytr_pu)*YI  */
						}
						break;
				}

				for (j=0; j<3; j++){
					indexFromTo[j]   = 3*ext2int_busMap[trafo_fbus]-3+j;
					indexFromTo[j+3] = 3*ext2int_busMap[trafo_tbus]-3+j;
					equivBus[indexFromTo[j]].isConnBus   = 1;
					equivBus[indexFromTo[j+3]].isConnBus = 1;
					equivBus[indexFromTo[j]].gl   += basePower*(Ypp_real[3*j] + Ypp_real[3*j+1] + \
						Ypp_real[3*j+2] + Yps_real[3*j] + Yps_real[3*j+1] + Yps_real[3*j+2]);
					equivBus[indexFromTo[j]].bl   += basePower*(Ypp_imag[3*j] + Ypp_imag[3*j+1] + \
						Ypp_imag[3*j+2] + Yps_imag[3*j] + Yps_imag[3*j+1] + Yps_imag[3*j+2]);
					equivBus[indexFromTo[j+3]].gl += basePower*(Ysp_real[3*j] + Ysp_real[3*j+1] + \
						Ysp_real[3*j+2] + Yss_real[3*j] + Yss_real[3*j+1] + Yss_real[3*j+2]);
					equivBus[indexFromTo[j+3]].bl += basePower*(Ysp_imag[3*j] + Ysp_imag[3*j+1] + \
						Ysp_imag[3*j+2] + Yss_imag[3*j] + Yss_imag[3*j+1] + Yss_imag[3*j+2]);
				}

				m = 0;
				for (j=0; j<5; j++){
					for (k=j+1; k<6; k++){
						equivBranch[15*j5+m].internal_i = indexFromTo[j];
						equivBranch[15*j5+m].internal_j = indexFromTo[k];
						m++;
					}
				}

				lineCond_val[0]  = - Ypp_real[1];		lineSusc_val[0]  = - Ypp_imag[1];
				lineCond_val[1]  = - Ypp_real[2];		lineSusc_val[1]  = - Ypp_imag[2];
				lineCond_val[2]  = - Yps_real[0];		lineSusc_val[2]  = - Yps_imag[0];
				lineCond_val[3]  = - Yps_real[1];		lineSusc_val[3]  = - Yps_imag[1];
				lineCond_val[4]  = - Yps_real[2];		lineSusc_val[4]  = - Yps_imag[2];
				lineCond_val[5]  = - Ypp_real[5];		lineSusc_val[5]  = - Ypp_imag[5];
				lineCond_val[6]  = - Yps_real[3];		lineSusc_val[6]  = - Yps_imag[3];
				lineCond_val[7]  = - Yps_real[4];		lineSusc_val[7]  = - Yps_imag[4];
				lineCond_val[8]  = - Yps_real[5];		lineSusc_val[8]  = - Yps_imag[5];
				lineCond_val[9]  = - Yps_real[6];		lineSusc_val[9]  = - Yps_imag[6];
				lineCond_val[10] = - Yps_real[7];		lineSusc_val[10] = - Yps_imag[7];
				lineCond_val[11] = - Yps_real[8];		lineSusc_val[11] = - Yps_imag[8];
				lineCond_val[12] = - Yss_real[1];		lineSusc_val[12] = - Yss_imag[1];
				lineCond_val[13] = - Yss_real[2];		lineSusc_val[13] = - Yss_imag[2];
				lineCond_val[14] = - Yss_real[5];		lineSusc_val[14] = - Yss_imag[5];

				for (j=0; j<15; j++){
					equivBranch[15*j5+j].status = 1;
					equivBranch[15*j5+j].yff[0] = equivBranch[15*j5+j].ytt[0] = lineCond_val[j];
					equivBranch[15*j5+j].yff[1] = equivBranch[15*j5+j].ytt[1]	= lineSusc_val[j];
					equivBranch[15*j5+j].yft[0] = equivBranch[15*j5+j].ytf[0] = -lineCond_val[j];
					equivBranch[15*j5+j].yft[1] = equivBranch[15*j5+j].ytf[1] = -lineSusc_val[j];
				}

			} /* Trafo_status conditional */
			j5++;
		} /* End of transformer-branch data*/

	} /* End of for (file reading) */
	fclose(fp);

	/**************************************************** Final arrangements ****************************************************************
		- Phases not utilized from each bus is treated as an isolated bus
		- Generator and load data structures are reorder according to bus numbers
	****************************************************************************************************************************************/

	for (i=0; i<pf->nbranch; i++){
		if (!(equivBranch[i].yff[0] || equivBranch[i].yff[1] || equivBranch[i].ytt[0] || \
			equivBranch[i].ytt[1] || equivBranch[i].yft[0] || equivBranch[i].yft[1] || \
			equivBranch[i].ytf[0] || equivBranch[i].ytf[1])){
			equivBranch[i].status = 0;
		}
	}

	for (i=0; i<pf->nbus; i++){
		if (!equivBus[i].isConnBus){
			equivBus[i].ide = ISOLATED_BUS;
		}
	}

	j1=0; j2=0;
	ierr = PetscCalloc1(pf->ngen,&equivGen_ord);CHKERRQ(ierr);
	ierr = PetscCalloc1(pf->nwye,&equivWyeLoad_ord);CHKERRQ(ierr);
	for (i = 0; i < pf->nbus; i++){
		if(pf->bus[i].ngen == 1){
			ierr = PetscMemcpy(&equivGen_ord[j1++],&pf->gen[pf->bus[i].gidx[0]],sizeof(struct _p_GEN));CHKERRQ(ierr);
		}
		if(pf->bus[i].nload == 1){
			ierr = PetscMemcpy(&equivWyeLoad_ord[j2++],&pf->wye[pf->bus[i].lidx[0]],sizeof(struct _p_WYELOAD));CHKERRQ(ierr);
		}
	}
	ierr = PetscFree(pf->gen);CHKERRQ(ierr);
	ierr = PetscFree(pf->wye);CHKERRQ(ierr);
	pf->gen = equivGen_ord;
	pf->wye = equivWyeLoad_ord;

	ierr = PetscFree(ext2int_busMap);CHKERRQ(ierr);

	/* End PETSc function */
	PetscFunctionReturn(0);

} /* End of function DPFReadMatPowerData */

#if 0

#endif
