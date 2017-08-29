static char help[] = "Run this program:\n\
                      mpiexec -n 1 ./dpflow\n\
                      mpiexec -n 1 ./dpflow -pfdata <filename>\n";
                      
/* T
   Concepts: DMNetwork
   Concepts: PETSc SNES solver
*/

#include "pf.h"

PetscMPIInt rank;

#undef __FUNCT__
#define __FUNCT__ "GetListofEdges"
PetscErrorCode GetListofEdges(PetscInt nbranches, EDGEDATA branch,int edges[])
{
  PetscInt       i, fbus,tbus;

  PetscFunctionBegin;
  for (i=0; i < nbranches; i++) {
    fbus = branch[i].internal_i;
    tbus = branch[i].internal_j;
    edges[2*i]   = fbus;
    edges[2*i+1] = tbus;
  }
  PetscFunctionReturn(0);
}

typedef struct{
  PetscScalar  Sbase;
}UserCtx;

#undef __FUNCT__
#define __FUNCT__ "FormFunction"
PetscErrorCode FormFunction(SNES snes,Vec X, Vec F,void *appctx)
{
  PetscErrorCode ierr;
  DM             networkdm;
  UserCtx       *User=(UserCtx*)appctx;
  Vec           localX,localF;
  PetscInt      e;
  PetscInt      v,v3,vStart,vEnd,vfrom,vto;
  const PetscScalar *xarr;
  PetscScalar   *farr;
  PetscInt      offsetfrom,offsetto,offset;
  DMNetworkComponentGenericDataType *arr;
  PetscScalar   Vm_A, Vm_B, Vm_C, Vm_AB, Vm_BC, Vm_CA;
  PetscScalar   Va_A, Va_B, Va_C, Va_AB, Va_BC, Va_CA;
  PetscScalar   P_AB,  P_BC,  P_CA,  Q_AB,  Q_BC,  Q_CA;
  PetscScalar   IP_AB, IP_BC, IP_CA, IQ_AB, IQ_BC, IQ_CA;
  PetscInt      idxA, idxB, idxC;

  PetscFunctionBegin;
  ierr = SNESGetDM(snes,&networkdm);CHKERRQ(ierr);
  ierr = DMGetLocalVector(networkdm,&localX);CHKERRQ(ierr);
  ierr = DMGetLocalVector(networkdm,&localF);CHKERRQ(ierr);
  ierr = VecSet(F,0.0);CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(networkdm,F,INSERT_VALUES,localF);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(networkdm,F,INSERT_VALUES,localF);CHKERRQ(ierr);

  ierr = VecGetArrayRead(localX,&xarr);CHKERRQ(ierr);
  ierr = VecGetArray(localF,&farr);CHKERRQ(ierr);

  ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr);
  ierr = DMNetworkGetComponentDataArray(networkdm,&arr);CHKERRQ(ierr);

  for (v3=0; v3<(vEnd-vStart)/3; v3++){
    P_AB=0.0,  P_BC=0.0,  P_CA=0.0,  Q_AB=0.0,  Q_BC=0.0,  Q_CA=0.0;
    IP_AB=0.0, IP_BC=0.0, IP_CA=0.0, IQ_AB=0.0, IQ_BC=0.0, IQ_CA=0.0;
    for (v=vStart+3*v3; v<vStart+3*v3+3; v++){
      PetscInt    i,j,offsetd,key;
      PetscScalar Vm;
      PetscScalar Sbase=User->Sbase;
      VERTEXDATA  bus=NULL;
      GEN         gen;
      WYELOAD     wye;
      DELTALOAD   delta;
      PetscInt    numComps;

      ierr = DMNetworkGetNumComponents(networkdm,v,&numComps);CHKERRQ(ierr);
      ierr = DMNetworkGetVariableOffset(networkdm,v,&offset);CHKERRQ(ierr);

      for (j = 0; j < numComps; j++) {
        ierr = DMNetworkGetComponentTypeOffset(networkdm,v,j,&key,&offsetd);CHKERRQ(ierr);
        if (key == 1) {
          PetscInt       nconnedges;
          const PetscInt *connedges;
  	      bus = (VERTEXDATA)(arr+offsetd);
        	/* Handle reference bus constrained dofs */
        	if (bus->ide == REF_BUS || bus->ide == ISOLATED_BUS) {
            farr[offset] = xarr[offset] - bus->va*PETSC_PI/180.0;
        	  farr[offset+1] = xarr[offset+1] - bus->vm;
            break;
        	}

      	  /* Shunt injections s*/
          Vm = xarr[offset+1];
      	  farr[offset] += Vm*Vm*bus->gl/Sbase;
      	  if(bus->ide != PV_BUS) farr[offset+1] += -Vm*Vm*bus->bl/Sbase;

  	      ierr = DMNetworkGetSupportingEdges(networkdm,v,&nconnedges,&connedges);CHKERRQ(ierr);
        	for (i=0; i < nconnedges; i++) {
        	  EDGEDATA       branch;
        	  PetscInt       keye;
            PetscScalar    Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
            const PetscInt *cone;
            PetscScalar    Vmf,Vmt,thetaf,thetat,thetaft,thetatf;

        	  e = connedges[i];
        	  ierr = DMNetworkGetComponentTypeOffset(networkdm,e,0,&keye,&offsetd);CHKERRQ(ierr);
        	  branch = (EDGEDATA)(arr+offsetd);
        	  if (!branch->status) continue;
        	  Gff = branch->yff[0];
        	  Bff = branch->yff[1];
        	  Gft = branch->yft[0];
        	  Bft = branch->yft[1];
        	  Gtf = branch->ytf[0];
        	  Btf = branch->ytf[1];
        	  Gtt = branch->ytt[0];
        	  Btt = branch->ytt[1];

        	  ierr = DMNetworkGetConnectedNodes(networkdm,e,&cone);CHKERRQ(ierr);
        	  vfrom = cone[0];
        	  vto   = cone[1];

        	  ierr = DMNetworkGetVariableOffset(networkdm,vfrom,&offsetfrom);CHKERRQ(ierr);
        	  ierr = DMNetworkGetVariableOffset(networkdm,vto,&offsetto);CHKERRQ(ierr);

        	  thetaf = xarr[offsetfrom];
        	  Vmf     = xarr[offsetfrom+1];
        	  thetat = xarr[offsetto];
        	  Vmt     = xarr[offsetto+1];
        	  thetaft = thetaf - thetat;
        	  thetatf = thetat - thetaf;

        	  if (vfrom == v) {
        	    farr[offsetfrom]   += Gff*Vmf*Vmf + Vmf*Vmt*(Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft));
        	    farr[offsetfrom+1] += -Bff*Vmf*Vmf + Vmf*Vmt*(-Bft*PetscCosScalar(thetaft) + Gft*PetscSinScalar(thetaft));
        	  } else {
        	    farr[offsetto]   += Gtt*Vmt*Vmt + Vmt*Vmf*(Gtf*PetscCosScalar(thetatf) + Btf*PetscSinScalar(thetatf));
        	    farr[offsetto+1] += -Btt*Vmt*Vmt + Vmt*Vmf*(-Btf*PetscCosScalar(thetatf) + Gtf*PetscSinScalar(thetatf));
        	  }
        	}
        } else if (key == 2) {
        	gen = (GEN)(arr+offsetd);
      	  farr[offset] += -gen->pg/Sbase;
      	  farr[offset+1] += -gen->qg/Sbase;
        } else if (key == 3) {
        	wye = (WYELOAD)(arr+offsetd);
          farr[offset] +=   (wye->pl + xarr[offset+1]*wye->ip)/Sbase;
          farr[offset+1] += (wye->ql + xarr[offset+1]*wye->iq)/Sbase;;
        } else if (key == 4) {
          delta = (DELTALOAD)(arr+offsetd);
          if ((v - vStart)%3 == 0){ /*Phase A*/
            P_AB  = delta->pl/Sbase;
            Q_AB  = delta->ql/Sbase;
            IP_AB = delta->ip/Sbase;
            IQ_AB = delta->iq/Sbase;
            idxA  = offset;
            Va_A  = xarr[offset];
            Vm_A = xarr[offset+1];
          } else if ((v - vStart)%3 == 1){ /*Phase B*/
            P_BC  = delta->pl/Sbase;
            Q_BC  = delta->ql/Sbase;
            IP_BC = delta->ip/Sbase;
            IQ_BC = delta->iq/Sbase;
            idxB  = offset;
            Va_B  = xarr[offset];
            Vm_B  = xarr[offset+1];
          } else { /*Phase C*/
            P_CA  = delta->pl/Sbase;
            Q_CA  = delta->ql/Sbase;
            IP_CA = delta->ip/Sbase;
            IQ_CA = delta->iq/Sbase;
            idxC  = offset;
            Va_C  = xarr[offset];
            Vm_C  = xarr[offset+1];
          }
        } /*End of key == 4*/

      } /* End of components per phase per vertex loop */
      if (bus && bus->ide == PV_BUS) {
        farr[offset+1] = xarr[offset+1] - bus->vm;
      }
    } /* End of phases-per-vertex for loop */

    if(P_AB || Q_AB || P_BC || Q_BC || P_CA || Q_CA || IP_AB || IQ_AB || IP_BC || IQ_BC || IP_CA || IQ_CA){
      Vm_AB = PetscSqrtScalar(Vm_A*Vm_A + Vm_B*Vm_B - 2*Vm_A*Vm_B*PetscCosScalar(Va_A - Va_B));
      Vm_BC = PetscSqrtScalar(Vm_B*Vm_B + Vm_C*Vm_C - 2*Vm_B*Vm_C*PetscCosScalar(Va_B - Va_C));
      Vm_CA = PetscSqrtScalar(Vm_C*Vm_C + Vm_A*Vm_A - 2*Vm_C*Vm_A*PetscCosScalar(Va_C - Va_A));
      if((Vm_A*PetscCosScalar(Va_A) - Vm_B*PetscCosScalar(Va_B)) > 0){
        Va_AB = PetscAtanReal((Vm_A*PetscSinScalar(Va_A) - Vm_B*PetscSinScalar(Va_B))/(Vm_A*PetscCosScalar(Va_A) - Vm_B*PetscCosScalar(Va_B)));
      } else{
        Va_AB = PETSC_PI + PetscAtanReal((Vm_A*PetscSinScalar(Va_A) - Vm_B*PetscSinScalar(Va_B))/(Vm_A*PetscCosScalar(Va_A) - Vm_B*PetscCosScalar(Va_B)));
      }
      if((Vm_B*PetscCosScalar(Va_B) - Vm_C*PetscCosScalar(Va_C)) > 0){
        Va_BC = PetscAtanReal((Vm_B*PetscSinScalar(Va_B) - Vm_C*PetscSinScalar(Va_C))/(Vm_B*PetscCosScalar(Va_B) - Vm_C*PetscCosScalar(Va_C)));
      } else{
        Va_BC = PETSC_PI + PetscAtanReal((Vm_B*PetscSinScalar(Va_B) - Vm_C*PetscSinScalar(Va_C))/(Vm_B*PetscCosScalar(Va_B) - Vm_C*PetscCosScalar(Va_C)));
      }
      if((Vm_C*PetscCosScalar(Va_C) - Vm_A*PetscCosScalar(Va_A)) > 0){
        Va_CA = PetscAtanReal((Vm_C*PetscSinScalar(Va_C) - Vm_A*PetscSinScalar(Va_A))/(Vm_C*PetscCosScalar(Va_C) - Vm_A*PetscCosScalar(Va_A)));
      } else{
        Va_CA = PETSC_PI + PetscAtanReal((Vm_C*PetscSinScalar(Va_C) - Vm_A*PetscSinScalar(Va_A))/(Vm_C*PetscCosScalar(Va_C) - Vm_A*PetscCosScalar(Va_A)));
      }

      if(P_AB || Q_AB || P_BC || Q_BC || P_CA || Q_CA){
        farr[idxA]   += + (Vm_A/Vm_AB) * (P_AB*PetscCosScalar(Va_A - Va_AB) - Q_AB*PetscSinScalar(Va_A - Va_AB))\
                        - (Vm_A/Vm_CA) * (P_CA*PetscCosScalar(Va_A - Va_CA) - Q_CA*PetscSinScalar(Va_A - Va_CA)); /*PA*/
        farr[idxA+1] += + (Vm_A/Vm_AB) * (P_AB*PetscSinScalar(Va_A - Va_AB) + Q_AB*PetscCosScalar(Va_A - Va_AB))\
                        - (Vm_A/Vm_CA) * (P_CA*PetscSinScalar(Va_A - Va_CA) + Q_CA*PetscCosScalar(Va_A - Va_CA)); /*QA*/
        farr[idxB]   += + (Vm_B/Vm_BC) * (P_BC*PetscCosScalar(Va_B - Va_BC) - Q_BC*PetscSinScalar(Va_B - Va_BC))\
                        - (Vm_B/Vm_AB) * (P_AB*PetscCosScalar(Va_B - Va_AB) - Q_AB*PetscSinScalar(Va_B - Va_AB)); /*PB*/
        farr[idxB+1] += + (Vm_B/Vm_BC) * (P_BC*PetscSinScalar(Va_B - Va_BC) + Q_BC*PetscCosScalar(Va_B - Va_BC))\
                        - (Vm_B/Vm_AB) * (P_AB*PetscSinScalar(Va_B - Va_AB) + Q_AB*PetscCosScalar(Va_B - Va_AB)); /*QB*/
        farr[idxC]   += + (Vm_C/Vm_CA) * (P_CA*PetscCosScalar(Va_C - Va_CA) - Q_CA*PetscSinScalar(Va_C - Va_CA))\
                        - (Vm_C/Vm_BC) * (P_BC*PetscCosScalar(Va_C - Va_BC) - Q_BC*PetscSinScalar(Va_C - Va_BC)); /*PC*/
        farr[idxC+1] += + (Vm_C/Vm_CA) * (P_CA*PetscSinScalar(Va_C - Va_CA) + Q_CA*PetscCosScalar(Va_C - Va_CA))\
                        - (Vm_C/Vm_BC) * (P_BC*PetscSinScalar(Va_C - Va_BC) + Q_BC*PetscCosScalar(Va_C - Va_BC)); /*QC*/
      }
      if(IP_AB || IQ_AB || IP_BC || IQ_BC || IP_CA || IQ_CA){
        farr[idxA]   += + (Vm_A/PetscSqrtScalar(3)) * (IP_AB*PetscCosScalar(Va_A - Va_AB) - IQ_AB*PetscSinScalar(Va_A - Va_AB))\
                        - (Vm_A/PetscSqrtScalar(3)) * (IP_CA*PetscCosScalar(Va_A - Va_CA) - IQ_CA*PetscSinScalar(Va_A - Va_CA)); /*PA*/
        farr[idxA+1] += + (Vm_A/PetscSqrtScalar(3)) * (IP_AB*PetscSinScalar(Va_A - Va_AB) + IQ_AB*PetscCosScalar(Va_A - Va_AB))\
                        - (Vm_A/PetscSqrtScalar(3)) * (IP_CA*PetscSinScalar(Va_A - Va_CA) + IQ_CA*PetscCosScalar(Va_A - Va_CA)); /*QA*/
        farr[idxB]   += + (Vm_B/PetscSqrtScalar(3)) * (IP_BC*PetscCosScalar(Va_B - Va_BC) - IQ_BC*PetscSinScalar(Va_B - Va_BC))\
                        - (Vm_B/PetscSqrtScalar(3)) * (IP_AB*PetscCosScalar(Va_B - Va_AB) - IQ_AB*PetscSinScalar(Va_B - Va_AB)); /*PB*/
        farr[idxB+1] += + (Vm_B/PetscSqrtScalar(3)) * (IP_BC*PetscSinScalar(Va_B - Va_BC) + IQ_BC*PetscCosScalar(Va_B - Va_BC))\
                        - (Vm_B/PetscSqrtScalar(3)) * (IP_AB*PetscSinScalar(Va_B - Va_AB) + IQ_AB*PetscCosScalar(Va_B - Va_AB)); /*QB*/
        farr[idxC]   += + (Vm_C/PetscSqrtScalar(3)) * (IP_CA*PetscCosScalar(Va_C - Va_CA) - IQ_CA*PetscSinScalar(Va_C - Va_CA))\
                        - (Vm_C/PetscSqrtScalar(3)) * (IP_BC*PetscCosScalar(Va_C - Va_BC) - IQ_BC*PetscSinScalar(Va_C - Va_BC)); /*PC*/
        farr[idxC+1] += + (Vm_C/PetscSqrtScalar(3)) * (IP_CA*PetscSinScalar(Va_C - Va_CA) + IQ_CA*PetscCosScalar(Va_C - Va_CA))\
                        - (Vm_C/PetscSqrtScalar(3)) * (IP_BC*PetscSinScalar(Va_C - Va_BC) + IQ_BC*PetscCosScalar(Va_C - Va_BC)); /*QC*/
      }

    } /* End of objective function addition caused by delta loads*/

  } /* End of vertices for loop */
  ierr = VecRestoreArrayRead(localX,&xarr);CHKERRQ(ierr);
  ierr = VecRestoreArray(localF,&farr);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(networkdm,&localX);CHKERRQ(ierr);

  ierr = DMLocalToGlobalBegin(networkdm,localF,ADD_VALUES,F);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(networkdm,localF,ADD_VALUES,F);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(networkdm,&localF);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormJacobian"
PetscErrorCode FormJacobian(SNES snes,Vec X, Mat J,Mat Jpre,void *appctx)
{
  PetscErrorCode ierr;
  DM            networkdm;
  UserCtx       *User=(UserCtx*)appctx;
  Vec           localX;
  PetscInt      e,v,v3,vStart,vEnd,vfrom,vto;
  const PetscScalar *xarr;
  PetscInt      offsetfrom,offsetto,goffsetfrom,goffsetto;
  DMNetworkComponentGenericDataType *arr;
  PetscInt      row[6],col[6];
  PetscScalar   values[36];
  PetscScalar   Vm_A, Vm_B, Vm_C, Vm_AB, Vm_BC, Vm_CA;
  PetscScalar   Va_A, Va_B, Va_C, Va_AB, Va_BC, Va_CA;
  PetscScalar   P_AB,  P_BC,  P_CA,  Q_AB,  Q_BC,  Q_CA;
  PetscScalar   IP_AB, IP_BC, IP_CA, IQ_AB, IQ_BC, IQ_CA;
  PetscInt      idxA, idxB, idxC;

  PetscFunctionBegin;
  ierr = MatZeroEntries(J);CHKERRQ(ierr);

  ierr = SNESGetDM(snes,&networkdm);CHKERRQ(ierr);
  ierr = DMGetLocalVector(networkdm,&localX);CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  ierr = VecGetArrayRead(localX,&xarr);CHKERRQ(ierr);

  ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr);
  ierr = DMNetworkGetComponentDataArray(networkdm,&arr);CHKERRQ(ierr);

  for (v3=0; v3<(vEnd-vStart)/3; v3++){
    P_AB=0.0,  P_BC=0.0,  P_CA=0.0,  Q_AB=0.0,  Q_BC=0.0,  Q_CA=0.0;
    IP_AB=0.0, IP_BC=0.0, IP_CA=0.0, IQ_AB=0.0, IQ_BC=0.0, IQ_CA=0.0;
    for (v=vStart+3*v3; v<vStart+3*v3+3; v++){
      PetscInt    i,j,offsetd,key;
      PetscInt    offset,goffset;
      PetscScalar Sbase=User->Sbase;
      VERTEXDATA  bus;
      WYELOAD     wye;
      DELTALOAD   delta;
      PetscInt    numComps;

      ierr = DMNetworkGetNumComponents(networkdm,v,&numComps);CHKERRQ(ierr);
      ierr = DMNetworkGetVariableOffset(networkdm,v,&offset);CHKERRQ(ierr);
      ierr = DMNetworkGetVariableGlobalOffset(networkdm,v,&goffset);CHKERRQ(ierr);

      for (j = 0; j < numComps; j++) {
        ierr = DMNetworkGetComponentTypeOffset(networkdm,v,j,&key,&offsetd);CHKERRQ(ierr);
        if (key == 1) {
          PetscInt       nconnedges;
        	const PetscInt *connedges;

        	bus = (VERTEXDATA)(arr+offsetd);

      	  /* Handle reference bus constrained dofs */
      	  if (bus->ide == REF_BUS || bus->ide == ISOLATED_BUS) {
      	    row[0] = goffset; row[1] = goffset+1;
      	    col[0] = goffset; col[1] = goffset+1;
      	    values[0] = 1.0; values[1] = 0.0; values[2] = 0.0; values[3] = 1.0;
      	    ierr = MatSetValues(J,2,row,2,col,values,ADD_VALUES);CHKERRQ(ierr);
      	    break;
      	  }

      	  /* Shunt injections */
          row[0] = goffset; row[1] = goffset+1;
          col[0] = goffset; col[1] = goffset+1;
          values[0] = values[1] = values[2] = values[3] = 0.0;
          if (bus->ide != PV_BUS) {
            values[1] = 2.0*xarr[offset+1]*bus->gl/Sbase;
            values[3] = -2.0*xarr[offset+1]*bus->bl/Sbase;
          }
          ierr = MatSetValues(J,2,row,2,col,values,ADD_VALUES);CHKERRQ(ierr);

        	ierr = DMNetworkGetSupportingEdges(networkdm,v,&nconnedges,&connedges);CHKERRQ(ierr);
        	for (i=0; i < nconnedges; i++) {
        	  EDGEDATA       branch;
        	  VERTEXDATA     busf,bust;
        	  PetscInt       offsetfd,offsettd,keyf,keyt;
            PetscScalar    Gff,Bff,Gft,Bft,Gtf,Btf,Gtt,Btt;
            const PetscInt *cone;
            PetscScalar    Vmf,Vmt,thetaf,thetat,thetaft,thetatf;

        	  e = connedges[i];
        	  ierr = DMNetworkGetComponentTypeOffset(networkdm,e,0,&key,&offsetd);CHKERRQ(ierr);
        	  branch = (EDGEDATA)(arr+offsetd);
        	  if (!branch->status) continue;

        	  Gff = branch->yff[0];
        	  Bff = branch->yff[1];
        	  Gft = branch->yft[0];
        	  Bft = branch->yft[1];
        	  Gtf = branch->ytf[0];
        	  Btf = branch->ytf[1];
        	  Gtt = branch->ytt[0];
        	  Btt = branch->ytt[1];

        	  ierr = DMNetworkGetConnectedNodes(networkdm,e,&cone);CHKERRQ(ierr);
        	  vfrom = cone[0];
        	  vto   = cone[1];

        	  ierr = DMNetworkGetVariableOffset(networkdm,vfrom,&offsetfrom);CHKERRQ(ierr);
        	  ierr = DMNetworkGetVariableOffset(networkdm,vto,&offsetto);CHKERRQ(ierr);
        	  ierr = DMNetworkGetVariableGlobalOffset(networkdm,vfrom,&goffsetfrom);CHKERRQ(ierr);
        	  ierr = DMNetworkGetVariableGlobalOffset(networkdm,vto,&goffsetto);CHKERRQ(ierr);

        	  if (goffsetto < 0) goffsetto = -goffsetto - 1;

        	  thetaf = xarr[offsetfrom];
        	  Vmf     = xarr[offsetfrom+1];
        	  thetat = xarr[offsetto];
        	  Vmt     = xarr[offsetto+1];
        	  thetaft = thetaf - thetat;
        	  thetatf = thetat - thetaf;

        	  ierr = DMNetworkGetComponentTypeOffset(networkdm,vfrom,0,&keyf,&offsetfd);CHKERRQ(ierr);
        	  ierr = DMNetworkGetComponentTypeOffset(networkdm,vto,0,&keyt,&offsettd);CHKERRQ(ierr);
        	  busf = (VERTEXDATA)(arr+offsetfd);
        	  bust = (VERTEXDATA)(arr+offsettd);

        	  if (vfrom == v) {
        	    if (busf->ide != REF_BUS) {
        	      /*    farr[offsetfrom]   += Gff*Vmf*Vmf + Vmf*Vmt*(Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft));  */
        	      row[0]  = goffsetfrom;
        	      col[0]  = goffsetfrom; col[1] = goffsetfrom+1; col[2] = goffsetto; col[3] = goffsetto+1;
        	      values[0] =  Vmf*Vmt*(Gft*-PetscSinScalar(thetaft) + Bft*PetscCosScalar(thetaft)); /* df_dthetaf */
        	      values[1] =  2.0*Gff*Vmf + Vmt*(Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft)); /* df_dVmf */
        	      values[2] =  Vmf*Vmt*(Gft*PetscSinScalar(thetaft) + Bft*-PetscCosScalar(thetaft)); /* df_dthetat */
        	      values[3] =  Vmf*(Gft*PetscCosScalar(thetaft) + Bft*PetscSinScalar(thetaft)); /* df_dVmt */

        	      ierr = MatSetValues(J,1,row,4,col,values,ADD_VALUES);CHKERRQ(ierr);
        	    }
        	    if (busf->ide != PV_BUS && busf->ide != REF_BUS) {
        	      row[0] = goffsetfrom+1;
        	      col[0]  = goffsetfrom; col[1] = goffsetfrom+1; col[2] = goffsetto; col[3] = goffsetto+1;
        	      /*    farr[offsetfrom+1] += -Bff*Vmf*Vmf + Vmf*Vmt*(-Bft*PetscCosScalar(thetaft) + Gft*PetscSinScalar(thetaft)); */
        	      values[0] =  Vmf*Vmt*(Bft*PetscSinScalar(thetaft) + Gft*PetscCosScalar(thetaft));
        	      values[1] =  -2.0*Bff*Vmf + Vmt*(-Bft*PetscCosScalar(thetaft) + Gft*PetscSinScalar(thetaft));
        	      values[2] =  Vmf*Vmt*(-Bft*PetscSinScalar(thetaft) + Gft*-PetscCosScalar(thetaft));
        	      values[3] =  Vmf*(-Bft*PetscCosScalar(thetaft) + Gft*PetscSinScalar(thetaft));

        	      ierr = MatSetValues(J,1,row,4,col,values,ADD_VALUES);CHKERRQ(ierr);
        	    }
        	  } else {
        	    if (bust->ide != REF_BUS) {
        	      row[0] = goffsetto;
        	      col[0] = goffsetto; col[1] = goffsetto+1; col[2] = goffsetfrom; col[3] = goffsetfrom+1;
        	      /*    farr[offsetto]   += Gtt*Vmt*Vmt + Vmt*Vmf*(Gtf*PetscCosScalar(thetatf) + Btf*PetscSinScalar(thetatf)); */
        	      values[0] =  Vmt*Vmf*(Gtf*-PetscSinScalar(thetatf) + Btf*PetscCosScalar(thetaft)); /* df_dthetat */
        	      values[1] =  2.0*Gtt*Vmt + Vmf*(Gtf*PetscCosScalar(thetatf) + Btf*PetscSinScalar(thetatf)); /* df_dVmt */
        	      values[2] =  Vmt*Vmf*(Gtf*PetscSinScalar(thetatf) + Btf*-PetscCosScalar(thetatf)); /* df_dthetaf */
        	      values[3] =  Vmt*(Gtf*PetscCosScalar(thetatf) + Btf*PetscSinScalar(thetatf)); /* df_dVmf */

        	      ierr = MatSetValues(J,1,row,4,col,values,ADD_VALUES);CHKERRQ(ierr);
        	    }
        	    if (bust->ide != PV_BUS && bust->ide != REF_BUS) {
        	      row[0] = goffsetto+1;
        	      col[0] = goffsetto; col[1] = goffsetto+1; col[2] = goffsetfrom; col[3] = goffsetfrom+1;
        	      /*    farr[offsetto+1] += -Btt*Vmt*Vmt + Vmt*Vmf*(-Btf*PetscCosScalar(thetatf) + Gtf*PetscSinScalar(thetatf)); */
        	      values[0] =  Vmt*Vmf*(Btf*PetscSinScalar(thetatf) + Gtf*PetscCosScalar(thetatf));
        	      values[1] =  -2.0*Btt*Vmt + Vmf*(-Btf*PetscCosScalar(thetatf) + Gtf*PetscSinScalar(thetatf));
        	      values[2] =  Vmt*Vmf*(-Btf*PetscSinScalar(thetatf) + Gtf*-PetscCosScalar(thetatf));
        	      values[3] =  Vmt*(-Btf*PetscCosScalar(thetatf) + Gtf*PetscSinScalar(thetatf));

        	      ierr = MatSetValues(J,1,row,4,col,values,ADD_VALUES);CHKERRQ(ierr);
        	    }
        	  }
        	}
        } else if (key == 3) {
      	  wye = (WYELOAD)(arr+offsetd);
          if (wye->ip/Sbase!=0.0 || wye->iq/Sbase!=0.0){
            row[0] = goffset; row[1] = goffset+1;
            col[0] = goffset+1;
            values[0] = wye->ip/Sbase;
            values[1] = wye->iq/Sbase;
            ierr = MatSetValues(J,2,row,1,col,values,ADD_VALUES);CHKERRQ(ierr);
          }
        } else if (key == 4) {
          delta = (DELTALOAD)(arr+offsetd);
          if ((v - vStart)%3 == 0){ /*Phase A*/
            P_AB  = delta->pl/Sbase;
            Q_AB  = delta->ql/Sbase;
            IP_AB = delta->ip/Sbase;
            IQ_AB = delta->iq/Sbase;
            idxA  = goffset;
            Va_A  = xarr[offset];
            Vm_A  = xarr[offset+1];
          } else if ((v - vStart)%3 == 1){ /*Phase B*/
            P_BC  = delta->pl/Sbase;
            Q_BC  = delta->ql/Sbase;
            IP_BC = delta->ip/Sbase;
            IQ_BC = delta->iq/Sbase;
            idxB  = goffset;
            Va_B  = xarr[offset];
            Vm_B  = xarr[offset+1];
          } else { /*Phase C*/
            P_CA  = delta->pl/Sbase;
            Q_CA  = delta->ql/Sbase;
            IP_CA = delta->ip/Sbase;
            IQ_CA = delta->iq/Sbase;
            idxC  = goffset;
            Va_C  = xarr[offset];
            Vm_C  = xarr[offset+1];
          }
        } /*End of key == 4*/

      } /*End of components per vertex loop*/
      if (bus && bus->ide == PV_BUS) {
        row[0] = goffset+1; col[0] = goffset+1;
        values[0]  = 1.0;
        ierr = MatSetValues(J,1,row,1,col,values,ADD_VALUES);CHKERRQ(ierr);
      }
    } /* End of phase-per-vertex for loop */

    if (P_AB || Q_AB || P_BC || Q_BC || P_CA || Q_CA || IP_AB || IQ_AB || IP_BC || IQ_BC || IP_CA || IQ_CA){
      Vm_AB = PetscSqrtScalar(Vm_A*Vm_A + Vm_B*Vm_B - 2*Vm_A*Vm_B*PetscCosScalar(Va_A - Va_B));
      Vm_BC = PetscSqrtScalar(Vm_B*Vm_B + Vm_C*Vm_C - 2*Vm_B*Vm_C*PetscCosScalar(Va_B - Va_C));
      Vm_CA = PetscSqrtScalar(Vm_C*Vm_C + Vm_A*Vm_A - 2*Vm_C*Vm_A*PetscCosScalar(Va_C - Va_A));
      if((Vm_A*PetscCosScalar(Va_A) - Vm_B*PetscCosScalar(Va_B)) > 0){
        Va_AB = PetscAtanReal((Vm_A*PetscSinScalar(Va_A) - Vm_B*PetscSinScalar(Va_B))/(Vm_A*PetscCosScalar(Va_A) - Vm_B*PetscCosScalar(Va_B)));
      } else{
        Va_AB = PETSC_PI + PetscAtanReal((Vm_A*PetscSinScalar(Va_A) - Vm_B*PetscSinScalar(Va_B))/(Vm_A*PetscCosScalar(Va_A) - Vm_B*PetscCosScalar(Va_B)));
      }
      if((Vm_B*PetscCosScalar(Va_B) - Vm_C*PetscCosScalar(Va_C)) > 0){
        Va_BC = PetscAtanReal((Vm_B*PetscSinScalar(Va_B) - Vm_C*PetscSinScalar(Va_C))/(Vm_B*PetscCosScalar(Va_B) - Vm_C*PetscCosScalar(Va_C)));
      } else{
        Va_BC = PETSC_PI + PetscAtanReal((Vm_B*PetscSinScalar(Va_B) - Vm_C*PetscSinScalar(Va_C))/(Vm_B*PetscCosScalar(Va_B) - Vm_C*PetscCosScalar(Va_C)));
      }
      if((Vm_C*PetscCosScalar(Va_C) - Vm_A*PetscCosScalar(Va_A)) > 0){
        Va_CA = PetscAtanReal((Vm_C*PetscSinScalar(Va_C) - Vm_A*PetscSinScalar(Va_A))/(Vm_C*PetscCosScalar(Va_C) - Vm_A*PetscCosScalar(Va_A)));
      } else{
        Va_CA = PETSC_PI + PetscAtanReal((Vm_C*PetscSinScalar(Va_C) - Vm_A*PetscSinScalar(Va_A))/(Vm_C*PetscCosScalar(Va_C) - Vm_A*PetscCosScalar(Va_A)));
      }

      row[0] = col[0] = idxA;
      row[1] = col[1] = idxA+1;
      row[2] = col[2] = idxB;
      row[3] = col[3] = idxB+1;
      row[4] = col[4] = idxC;
      row[5] = col[5] = idxC+1;

      if (P_AB || Q_AB || P_BC || Q_BC || P_CA || Q_CA){
        /* dP/dVang */
        values[0]  = + (Vm_A*Vm_B/(Vm_AB*Vm_AB))*(+ P_AB*PetscSinScalar(Va_A+Va_B-2*Va_AB) + Q_AB*PetscCosScalar(Va_A+Va_B-2*Va_AB)) \
                     + (Vm_A*Vm_C/(Vm_CA*Vm_CA))*(+ P_CA*PetscSinScalar(Va_A+Va_C-2*Va_CA) + Q_CA*PetscCosScalar(Va_A+Va_C-2*Va_CA)); /*dP_A/dVang_A*/
        values[2]  = - (Vm_A*Vm_B/(Vm_AB*Vm_AB))*(+ P_AB*PetscSinScalar(Va_A+Va_B-2*Va_AB) + Q_AB*PetscCosScalar(Va_A+Va_B-2*Va_AB)); /*dP_A/dVang_B*/
        values[4]  = - (Vm_A*Vm_C/(Vm_CA*Vm_CA))*(+ P_CA*PetscSinScalar(Va_A+Va_C-2*Va_CA) + Q_CA*PetscCosScalar(Va_A+Va_C-2*Va_CA)); /*dP_A/dVang_C*/

        values[12] = - (Vm_A*Vm_B/(Vm_AB*Vm_AB))*(+ P_AB*PetscSinScalar(Va_A+Va_B-2*Va_AB) + Q_AB*PetscCosScalar(Va_A+Va_B-2*Va_AB)); /*dP_B/dVang_A*/
        values[14] = + (Vm_B*Vm_C/(Vm_BC*Vm_BC))*(+ P_BC*PetscSinScalar(Va_B+Va_C-2*Va_BC) + Q_BC*PetscCosScalar(Va_B+Va_C-2*Va_BC)) \
                     + (Vm_A*Vm_B/(Vm_AB*Vm_AB))*(+ P_AB*PetscSinScalar(Va_A+Va_B-2*Va_AB) + Q_AB*PetscCosScalar(Va_A+Va_B-2*Va_AB)); /*dP_B/dVang_B*/
        values[16] = - (Vm_B*Vm_C/(Vm_BC*Vm_BC))*(+ P_BC*PetscSinScalar(Va_B+Va_C-2*Va_BC) + Q_BC*PetscCosScalar(Va_B+Va_C-2*Va_BC)); /*dP_B/dVang_C*/

        values[24] = - (Vm_A*Vm_C/(Vm_CA*Vm_CA))*(+ P_CA*PetscSinScalar(Va_A+Va_C-2*Va_CA) + Q_CA*PetscCosScalar(Va_A+Va_C-2*Va_CA)); /*dP_C/dVang_A*/
        values[26] = - (Vm_B*Vm_C/(Vm_BC*Vm_BC))*(+ P_BC*PetscSinScalar(Va_B+Va_C-2*Va_BC) + Q_BC*PetscCosScalar(Va_B+Va_C-2*Va_BC)); /*dP_C/dVang_B*/
        values[28] = + (Vm_A*Vm_C/(Vm_CA*Vm_CA))*(+ P_CA*PetscSinScalar(Va_A+Va_C-2*Va_CA) + Q_CA*PetscCosScalar(Va_A+Va_C-2*Va_CA)) \
                     + (Vm_B*Vm_C/(Vm_BC*Vm_BC))*(+ P_BC*PetscSinScalar(Va_B+Va_C-2*Va_BC) + Q_BC*PetscCosScalar(Va_B+Va_C-2*Va_BC)); /*dP_C/dVang_C*/

        /* dQ/dVang */
        values[6]  = + (Vm_A*Vm_B/(Vm_AB*Vm_AB))*(- P_AB*PetscCosScalar(Va_A+Va_B-2*Va_AB) + Q_AB*PetscSinScalar(Va_A+Va_B-2*Va_AB)) \
                     + (Vm_A*Vm_C/(Vm_CA*Vm_CA))*(- P_CA*PetscCosScalar(Va_A+Va_C-2*Va_CA) + Q_CA*PetscSinScalar(Va_A+Va_C-2*Va_CA)); /*dQ_A/dVang_A*/
        values[8]  = - (Vm_A*Vm_B/(Vm_AB*Vm_AB))*(- P_AB*PetscCosScalar(Va_A+Va_B-2*Va_AB) + Q_AB*PetscSinScalar(Va_A+Va_B-2*Va_AB)); /*dQ_A/dVang_B*/
        values[10] = - (Vm_A*Vm_C/(Vm_CA*Vm_CA))*(- P_CA*PetscCosScalar(Va_A+Va_C-2*Va_CA) + Q_CA*PetscSinScalar(Va_A+Va_C-2*Va_CA)); /*dQ_A/dVang_C*/

        values[18] = - (Vm_A*Vm_B/(Vm_AB*Vm_AB))*(- P_AB*PetscCosScalar(Va_A+Va_B-2*Va_AB) + Q_AB*PetscSinScalar(Va_A+Va_B-2*Va_AB)); /*dQ_B/dVang_A*/
        values[20] = + (Vm_B*Vm_C/(Vm_BC*Vm_BC))*(- P_BC*PetscCosScalar(Va_B+Va_C-2*Va_BC) + Q_BC*PetscSinScalar(Va_B+Va_C-2*Va_BC)) \
                     + (Vm_A*Vm_B/(Vm_AB*Vm_AB))*(- P_AB*PetscCosScalar(Va_A+Va_B-2*Va_AB) + Q_AB*PetscSinScalar(Va_A+Va_B-2*Va_AB)); /*dQ_B/dVang_B*/
        values[22] = - (Vm_B*Vm_C/(Vm_BC*Vm_BC))*(- P_BC*PetscCosScalar(Va_B+Va_C-2*Va_BC) + Q_BC*PetscSinScalar(Va_B+Va_C-2*Va_BC)); /*dQ_B/dVang_C*/

        values[30] = - (Vm_A*Vm_C/(Vm_CA*Vm_CA))*(- P_CA*PetscCosScalar(Va_A+Va_C-2*Va_CA) + Q_CA*PetscSinScalar(Va_A+Va_C-2*Va_CA)); /*dQ_C/dVang_A*/
        values[32] = - (Vm_B*Vm_C/(Vm_BC*Vm_BC))*(- P_BC*PetscCosScalar(Va_B+Va_C-2*Va_BC) + Q_BC*PetscSinScalar(Va_B+Va_C-2*Va_BC)); /*dQ_C/dVang_B*/
        values[34] = + (Vm_A*Vm_C/(Vm_CA*Vm_CA))*(- P_CA*PetscCosScalar(Va_A+Va_C-2*Va_CA) + Q_CA*PetscSinScalar(Va_A+Va_C-2*Va_CA)) \
                     + (Vm_B*Vm_C/(Vm_BC*Vm_BC))*(- P_BC*PetscCosScalar(Va_B+Va_C-2*Va_BC) + Q_BC*PetscSinScalar(Va_B+Va_C-2*Va_BC)); /*dQ_C/dVang_C*/

        /* dP/dVmag */
        values[1]  = - (Vm_B/(Vm_AB*Vm_AB))*(P_AB*PetscCosScalar(Va_A+Va_B-2*Va_AB) - Q_AB*PetscSinScalar(Va_A+Va_B-2*Va_AB)) \
                     - (Vm_C/(Vm_CA*Vm_CA))*(P_CA*PetscCosScalar(Va_A+Va_C-2*Va_CA) - Q_CA*PetscSinScalar(Va_A+Va_C-2*Va_CA)); /*dP_A/dVmag_A*/
        values[3]  = + (Vm_A/(Vm_AB*Vm_AB))*(P_AB*PetscCosScalar(Va_A+Va_B-2*Va_AB) - Q_AB*PetscSinScalar(Va_A+Va_B-2*Va_AB)); /*dP_A/dVmag_B*/
        values[5]  = + (Vm_A/(Vm_CA*Vm_CA))*(P_CA*PetscCosScalar(Va_A+Va_C-2*Va_CA) - Q_CA*PetscSinScalar(Va_A+Va_C-2*Va_CA)); /*dP_A/dVmag_C*/

        values[13] = + (Vm_B/(Vm_AB*Vm_AB))*(P_AB*PetscCosScalar(Va_A+Va_B-2*Va_AB) - Q_AB*PetscSinScalar(Va_A+Va_B-2*Va_AB)); /*dP_B/dVmag_A*/
        values[15] = - (Vm_C/(Vm_BC*Vm_BC))*(P_BC*PetscCosScalar(Va_B+Va_C-2*Va_BC) - Q_BC*PetscSinScalar(Va_B+Va_C-2*Va_BC)) \
                     - (Vm_A/(Vm_AB*Vm_AB))*(P_AB*PetscCosScalar(Va_A+Va_B-2*Va_AB) - Q_AB*PetscSinScalar(Va_A+Va_B-2*Va_AB)); /*dP_B/dVmag_B*/
        values[17] = + (Vm_B/(Vm_BC*Vm_BC))*(P_BC*PetscCosScalar(Va_B+Va_C-2*Va_BC) - Q_BC*PetscSinScalar(Va_B+Va_C-2*Va_BC)); /*dP_B/dVmag_C*/

        values[25] = + (Vm_C/(Vm_CA*Vm_CA))*(P_CA*PetscCosScalar(Va_A+Va_C-2*Va_CA) - Q_CA*PetscSinScalar(Va_A+Va_C-2*Va_CA)); /*dP_C/dVmag_A*/
        values[27] = + (Vm_C/(Vm_BC*Vm_BC))*(P_BC*PetscCosScalar(Va_B+Va_C-2*Va_BC) - Q_BC*PetscSinScalar(Va_B+Va_C-2*Va_BC)); /*dP_C/dVmag_B*/
        values[29] = - (Vm_A/(Vm_CA*Vm_CA))*(P_CA*PetscCosScalar(Va_A+Va_C-2*Va_CA) - Q_CA*PetscSinScalar(Va_A+Va_C-2*Va_CA)) \
                     - (Vm_B/(Vm_BC*Vm_BC))*(P_BC*PetscCosScalar(Va_B+Va_C-2*Va_BC) - Q_BC*PetscSinScalar(Va_B+Va_C-2*Va_BC)); /*dP_C/dVmag_C*/

        /* dQ/dVmag */
        values[7]  = - (Vm_B/(Vm_AB*Vm_AB))*(P_AB*PetscSinScalar(Va_A+Va_B-2*Va_AB) + Q_AB*PetscCosScalar(Va_A+Va_B-2*Va_AB)) \
                     - (Vm_C/(Vm_CA*Vm_CA))*(P_CA*PetscSinScalar(Va_A+Va_C-2*Va_CA) + Q_CA*PetscCosScalar(Va_A+Va_C-2*Va_CA)); /*dQ_A/dVmag_A*/
        values[9]  = + (Vm_A/(Vm_AB*Vm_AB))*(P_AB*PetscSinScalar(Va_A+Va_B-2*Va_AB) + Q_AB*PetscCosScalar(Va_A+Va_B-2*Va_AB)); /*dQ_A/dVmag_B*/
        values[11] = + (Vm_A/(Vm_CA*Vm_CA))*(P_CA*PetscSinScalar(Va_A+Va_C-2*Va_CA) + Q_CA*PetscCosScalar(Va_A+Va_C-2*Va_CA)); /*dQ_A/dVmag_C*/

        values[19] = + (Vm_B/(Vm_AB*Vm_AB))*(P_AB*PetscSinScalar(Va_A+Va_B-2*Va_AB) + Q_AB*PetscCosScalar(Va_A+Va_B-2*Va_AB)); /*dQ_B/dVmag_A*/
        values[21] = - (Vm_C/(Vm_BC*Vm_BC))*(P_BC*PetscSinScalar(Va_B+Va_C-2*Va_BC) + Q_BC*PetscCosScalar(Va_B+Va_C-2*Va_BC)) \
                     - (Vm_A/(Vm_AB*Vm_AB))*(P_AB*PetscSinScalar(Va_A+Va_B-2*Va_AB) + Q_AB*PetscCosScalar(Va_A+Va_B-2*Va_AB)); /*dQ_B/dVmag_B*/
        values[23] = + (Vm_B/(Vm_BC*Vm_BC))*(P_BC*PetscSinScalar(Va_B+Va_C-2*Va_BC) + Q_BC*PetscCosScalar(Va_B+Va_C-2*Va_BC)); /*dQ_B/dVmag_C*/

        values[31] = + (Vm_C/(Vm_CA*Vm_CA))*(P_CA*PetscSinScalar(Va_A+Va_C-2*Va_CA) + Q_CA*PetscCosScalar(Va_A+Va_C-2*Va_CA)); /*dQ_C/dVmag_A*/
        values[33] = + (Vm_C/(Vm_BC*Vm_BC))*(P_BC*PetscSinScalar(Va_B+Va_C-2*Va_BC) + Q_BC*PetscCosScalar(Va_B+Va_C-2*Va_BC)); /*dQ_C/dVmag_B*/
        values[35] = - (Vm_A/(Vm_CA*Vm_CA))*(P_CA*PetscSinScalar(Va_A+Va_C-2*Va_CA) + Q_CA*PetscCosScalar(Va_A+Va_C-2*Va_CA)) \
                     - (Vm_B/(Vm_BC*Vm_BC))*(P_BC*PetscSinScalar(Va_B+Va_C-2*Va_BC) + Q_BC*PetscCosScalar(Va_B+Va_C-2*Va_BC)); /*dQ_C/dVmag_C*/

        ierr = MatSetValues(J,6,row,6,col,values,ADD_VALUES);CHKERRQ(ierr);

      } /* End of conditional for constant-power delta loads */

      if (IP_AB || IQ_AB || IP_BC || IQ_BC || IP_CA || IQ_CA){
        /* dP/dVang */
        values[0]  = - (Vm_A*Vm_B/(PetscSqrtScalar(3.0)*Vm_AB))*(IP_AB*PetscCosScalar(Va_A+Va_B-2*Va_AB+PETSC_PI/2) - IQ_AB*PetscSinScalar(Va_A+Va_B-2*Va_AB+PETSC_PI/2))\
                     + (Vm_A/(PetscSqrtScalar(3.0)*Vm_AB*Vm_AB))*(Vm_A*Vm_B*PetscSinScalar(Va_A-Va_B))*(IP_AB*PetscCosScalar(Va_A-Va_AB) - IQ_AB*PetscSinScalar(Va_A-Va_AB))\
                     - (Vm_A*Vm_C/(PetscSqrtScalar(3.0)*Vm_CA))*(IP_CA*PetscCosScalar(Va_A+Va_C-2*Va_CA+PETSC_PI/2) - IQ_CA*PetscSinScalar(Va_A+Va_C-2*Va_CA+PETSC_PI/2))\
                     - (Vm_A/(PetscSqrtScalar(3.0)*Vm_CA*Vm_CA))*(Vm_A*Vm_C*PetscSinScalar(Va_A-Va_C))*(IP_CA*PetscCosScalar(Va_A-Va_CA) - IQ_CA*PetscSinScalar(Va_A-Va_CA)); /*dP_A/dVang_A*/
        values[2]  = + (Vm_B*Vm_A/(PetscSqrtScalar(3.0)*Vm_AB))*(IP_AB*PetscCosScalar(Va_A+Va_B-2*Va_AB+PETSC_PI/2) - IQ_AB*PetscSinScalar(Va_A+Va_B-2*Va_AB+PETSC_PI/2))\
                     + (Vm_A/(PetscSqrtScalar(3.0)*Vm_AB*Vm_AB))*(Vm_B*Vm_A*PetscSinScalar(Va_B-Va_A))*(IP_AB*PetscCosScalar(Va_A-Va_AB) - IQ_AB*PetscSinScalar(Va_A-Va_AB)); /*dP_A/dVang_B*/
        values[4]  = + (Vm_C*Vm_A/(PetscSqrtScalar(3.0)*Vm_CA))*(IP_CA*PetscCosScalar(Va_C+Va_A-2*Va_CA+PETSC_PI/2) - IQ_CA*PetscSinScalar(Va_C+Va_A-2*Va_CA+PETSC_PI/2))\
                     - (Vm_A/(PetscSqrtScalar(3.0)*Vm_CA*Vm_CA))*(Vm_C*Vm_A*PetscSinScalar(Va_C-Va_A))*(IP_CA*PetscCosScalar(Va_A-Va_CA) - IQ_CA*PetscSinScalar(Va_A-Va_CA)); /*dP_A/dVang_C*/

        values[12] = + (Vm_A*Vm_B/(PetscSqrtScalar(3.0)*Vm_AB))*(IP_AB*PetscCosScalar(Va_A+Va_B-2*Va_AB+PETSC_PI/2) - IQ_AB*PetscSinScalar(Va_A+Va_B-2*Va_AB+PETSC_PI/2))\
                     - (Vm_B/(PetscSqrtScalar(3.0)*Vm_AB*Vm_AB))*(Vm_A*Vm_B*PetscSinScalar(Va_A-Va_B))*(IP_AB*PetscCosScalar(Va_B-Va_AB) - IQ_AB*PetscSinScalar(Va_B-Va_AB)); /*dP_B/dVang_A*/
        values[14] = - (Vm_B*Vm_C/(PetscSqrtScalar(3.0)*Vm_BC))*(IP_BC*PetscCosScalar(Va_B+Va_C-2*Va_BC+PETSC_PI/2) - IQ_BC*PetscSinScalar(Va_B+Va_C-2*Va_BC+PETSC_PI/2))\
                     + (Vm_B/(PetscSqrtScalar(3.0)*Vm_BC*Vm_BC))*(Vm_B*Vm_C*PetscSinScalar(Va_B-Va_C))*(IP_BC*PetscCosScalar(Va_B-Va_BC) - IQ_BC*PetscSinScalar(Va_B-Va_BC))\
                     - (Vm_B*Vm_A/(PetscSqrtScalar(3.0)*Vm_AB))*(IP_AB*PetscCosScalar(Va_A+Va_B-2*Va_AB+PETSC_PI/2) - IQ_AB*PetscSinScalar(Va_A+Va_B-2*Va_AB+PETSC_PI/2))\
                     - (Vm_B/(PetscSqrtScalar(3.0)*Vm_AB*Vm_AB))*(Vm_B*Vm_A*PetscSinScalar(Va_B-Va_A))*(IP_AB*PetscCosScalar(Va_B-Va_AB) - IQ_AB*PetscSinScalar(Va_B-Va_AB)); /*dP_B/dVang_B*/
        values[16] = + (Vm_C*Vm_B/(PetscSqrtScalar(3.0)*Vm_BC))*(IP_BC*PetscCosScalar(Va_B+Va_C-2*Va_BC+PETSC_PI/2) - IQ_BC*PetscSinScalar(Va_B+Va_C-2*Va_BC+PETSC_PI/2))\
                     + (Vm_B/(PetscSqrtScalar(3.0)*Vm_BC*Vm_BC))*(Vm_C*Vm_B*PetscSinScalar(Va_C-Va_B))*(IP_BC*PetscCosScalar(Va_B-Va_BC) - IQ_BC*PetscSinScalar(Va_B-Va_BC)); /*dP_B/dVang_C*/

        values[24] = + (Vm_A*Vm_C/(PetscSqrtScalar(3.0)*Vm_CA))*(IP_CA*PetscCosScalar(Va_C+Va_A-2*Va_CA+PETSC_PI/2) - IQ_CA*PetscSinScalar(Va_C+Va_A-2*Va_CA+PETSC_PI/2))\
                     + (Vm_C/(PetscSqrtScalar(3.0)*Vm_CA*Vm_CA))*(Vm_A*Vm_C*PetscSinScalar(Va_A-Va_C))*(IP_CA*PetscCosScalar(Va_C-Va_CA) - IQ_CA*PetscSinScalar(Va_C-Va_CA)); /*dP_C/dVang_A*/
        values[26] = + (Vm_B*Vm_C/(PetscSqrtScalar(3.0)*Vm_BC))*(IP_BC*PetscCosScalar(Va_B+Va_C-2*Va_BC+PETSC_PI/2) - IQ_BC*PetscSinScalar(Va_B+Va_C-2*Va_BC+PETSC_PI/2))\
                     - (Vm_C/(PetscSqrtScalar(3.0)*Vm_BC*Vm_BC))*(Vm_B*Vm_C*PetscSinScalar(Va_B-Va_C))*(IP_BC*PetscCosScalar(Va_C-Va_BC) - IQ_BC*PetscSinScalar(Va_C-Va_BC)); /*dP_C/dVang_B*/
        values[28] = - (Vm_C*Vm_A/(PetscSqrtScalar(3.0)*Vm_CA))*(IP_CA*PetscCosScalar(Va_C+Va_A-2*Va_CA+PETSC_PI/2) - IQ_CA*PetscSinScalar(Va_C+Va_A-2*Va_CA+PETSC_PI/2))\
                     + (Vm_C/(PetscSqrtScalar(3.0)*Vm_CA*Vm_CA))*(Vm_C*Vm_A*PetscSinScalar(Va_C-Va_A))*(IP_CA*PetscCosScalar(Va_C-Va_CA) - IQ_CA*PetscSinScalar(Va_C-Va_CA))\
                     - (Vm_C*Vm_B/(PetscSqrtScalar(3.0)*Vm_BC))*(IP_BC*PetscCosScalar(Va_B+Va_C-2*Va_BC+PETSC_PI/2) - IQ_BC*PetscSinScalar(Va_B+Va_C-2*Va_BC+PETSC_PI/2))\
                     - (Vm_C/(PetscSqrtScalar(3.0)*Vm_BC*Vm_BC))*(Vm_C*Vm_B*PetscSinScalar(Va_C-Va_B))*(IP_BC*PetscCosScalar(Va_C-Va_BC) - IQ_BC*PetscSinScalar(Va_C-Va_BC)); /*dP_C/dVang_C*/

       /* dQ/dVang */
        values[6]  = - (Vm_A*Vm_B/(PetscSqrtScalar(3.0)*Vm_AB))*(IP_AB*PetscSinScalar(Va_A+Va_B-2*Va_AB+PETSC_PI/2) + IQ_AB*PetscCosScalar(Va_A+Va_B-2*Va_AB+PETSC_PI/2))\
                     + (Vm_A/(PetscSqrtScalar(3.0)*Vm_AB*Vm_AB))*(Vm_A*Vm_B*PetscSinScalar(Va_A-Va_B))*(IP_AB*PetscSinScalar(Va_A-Va_AB) + IQ_AB*PetscCosScalar(Va_A-Va_AB))\
                     - (Vm_A*Vm_C/(PetscSqrtScalar(3.0)*Vm_CA))*(IP_CA*PetscSinScalar(Va_A+Va_C-2*Va_CA+PETSC_PI/2) + IQ_CA*PetscCosScalar(Va_A+Va_C-2*Va_CA+PETSC_PI/2))\
                     - (Vm_A/(PetscSqrtScalar(3.0)*Vm_CA*Vm_CA))*(Vm_A*Vm_C*PetscSinScalar(Va_A-Va_C))*(IP_CA*PetscSinScalar(Va_A-Va_CA) + IQ_CA*PetscCosScalar(Va_A-Va_CA)); /*dQ_A/dVang_A*/
        values[8]  = + (Vm_B*Vm_A/(PetscSqrtScalar(3.0)*Vm_AB))*(IP_AB*PetscSinScalar(Va_A+Va_B-2*Va_AB+PETSC_PI/2) + IQ_AB*PetscCosScalar(Va_A+Va_B-2*Va_AB+PETSC_PI/2))\
                     + (Vm_A/(PetscSqrtScalar(3.0)*Vm_AB*Vm_AB))*(Vm_B*Vm_A*PetscSinScalar(Va_B-Va_A))*(IP_AB*PetscSinScalar(Va_A-Va_AB) + IQ_AB*PetscCosScalar(Va_A-Va_AB)); /*dQ_A/dVang_B*/
        values[10] = + (Vm_C*Vm_A/(PetscSqrtScalar(3.0)*Vm_CA))*(IP_CA*PetscSinScalar(Va_C+Va_A-2*Va_CA+PETSC_PI/2) + IQ_CA*PetscCosScalar(Va_C+Va_A-2*Va_CA+PETSC_PI/2))\
                     - (Vm_A/(PetscSqrtScalar(3.0)*Vm_CA*Vm_CA))*(Vm_C*Vm_A*PetscSinScalar(Va_C-Va_A))*(IP_CA*PetscSinScalar(Va_A - Va_CA) + IQ_CA*PetscCosScalar(Va_A - Va_CA)); /*dQ_A/dVang_C*/

        values[18] = + (Vm_A*Vm_B/(PetscSqrtScalar(3.0)*Vm_AB))*(IP_AB*PetscSinScalar(Va_A+Va_B-2*Va_AB+PETSC_PI/2) + IQ_AB*PetscCosScalar(Va_A+Va_B-2*Va_AB+PETSC_PI/2))\
                     - (Vm_B/(PetscSqrtScalar(3.0)*Vm_AB*Vm_AB))*(Vm_A*Vm_B*PetscSinScalar(Va_A-Va_B))*(IP_AB*PetscSinScalar(Va_B-Va_AB) + IQ_AB*PetscCosScalar(Va_B-Va_AB)); /*dQ_B/dVang_A*/
        values[20] = - (Vm_B*Vm_C/(PetscSqrtScalar(3.0)*Vm_BC))*(IP_BC*PetscSinScalar(Va_B+Va_C-2*Va_BC+PETSC_PI/2) + IQ_BC*PetscCosScalar(Va_B+Va_C-2*Va_BC+PETSC_PI/2))\
                     + (Vm_B/(PetscSqrtScalar(3.0)*Vm_BC*Vm_BC))*(Vm_B*Vm_C*PetscSinScalar(Va_B-Va_C))*(IP_BC*PetscSinScalar(Va_B-Va_BC) + IQ_BC*PetscCosScalar(Va_B-Va_BC))\
                     - (Vm_B*Vm_A/(PetscSqrtScalar(3.0)*Vm_AB))*(IP_AB*PetscSinScalar(Va_A+Va_B-2*Va_AB+PETSC_PI/2) + IQ_AB*PetscCosScalar(Va_A+Va_B-2*Va_AB+PETSC_PI/2))\
                     - (Vm_B/(PetscSqrtScalar(3.0)*Vm_AB*Vm_AB))*(Vm_B*Vm_A*PetscSinScalar(Va_B-Va_A))*(IP_AB*PetscSinScalar(Va_B-Va_AB) + IQ_AB*PetscCosScalar(Va_B-Va_AB)); /*dQ_B/dVang_B*/
        values[22] = + (Vm_C*Vm_B/(PetscSqrtScalar(3.0)*Vm_BC))*(IP_BC*PetscSinScalar(Va_B+Va_C-2*Va_BC+PETSC_PI/2) + IQ_BC*PetscCosScalar(Va_B+Va_C-2*Va_BC+PETSC_PI/2))\
                     + (Vm_B/(PetscSqrtScalar(3.0)*Vm_BC*Vm_BC))*(Vm_C*Vm_B*PetscSinScalar(Va_C-Va_B))*(IP_BC*PetscSinScalar(Va_B-Va_BC) + IQ_BC*PetscCosScalar(Va_B-Va_BC)); /*dQ_B/dVang_C*/

        values[30] = + (Vm_A*Vm_C/(PetscSqrtScalar(3.0)*Vm_CA))*(IP_CA*PetscSinScalar(Va_C+Va_A-2*Va_CA+PETSC_PI/2) + IQ_CA*PetscCosScalar(Va_C+Va_A-2*Va_CA+PETSC_PI/2))\
                     + (Vm_C/(PetscSqrtScalar(3.0)*Vm_CA*Vm_CA))*(Vm_A*Vm_C*PetscSinScalar(Va_A-Va_C))*(IP_CA*PetscSinScalar(Va_C-Va_CA) + IQ_CA*PetscCosScalar(Va_C-Va_CA)); /*dQ_C/dVang_A*/
        values[32] = + (Vm_B*Vm_C/(PetscSqrtScalar(3.0)*Vm_BC))*(IP_BC*PetscSinScalar(Va_B+Va_C-2*Va_BC+PETSC_PI/2) + IQ_BC*PetscCosScalar(Va_B+Va_C-2*Va_BC+PETSC_PI/2))\
                     - (Vm_C/(PetscSqrtScalar(3.0)*Vm_BC*Vm_BC))*(Vm_B*Vm_C*PetscSinScalar(Va_B-Va_C))*(IP_BC*PetscSinScalar(Va_C-Va_BC) + IQ_BC*PetscCosScalar(Va_C-Va_BC)); /*dQ_C/dVang_B*/
        values[34] = - (Vm_C*Vm_A/(PetscSqrtScalar(3.0)*Vm_CA))*(IP_CA*PetscSinScalar(Va_C+Va_A-2*Va_CA+PETSC_PI/2) + IQ_CA*PetscCosScalar(Va_C+Va_A-2*Va_CA+PETSC_PI/2))\
                     + (Vm_C/(PetscSqrtScalar(3.0)*Vm_CA*Vm_CA))*(Vm_C*Vm_A*PetscSinScalar(Va_C-Va_A))*(IP_CA*PetscSinScalar(Va_C-Va_CA) + IQ_CA*PetscCosScalar(Va_C-Va_CA))\
                     - (Vm_C*Vm_B/(PetscSqrtScalar(3.0)*Vm_BC))*(IP_BC*PetscSinScalar(Va_B+Va_C-2*Va_BC+PETSC_PI/2) + IQ_BC*PetscCosScalar(Va_B+Va_C-2*Va_BC+PETSC_PI/2))\
                     - (Vm_C/(PetscSqrtScalar(3.0)*Vm_BC*Vm_BC))*(Vm_C*Vm_B*PetscSinScalar(Va_C-Va_B))*(IP_BC*PetscSinScalar(Va_C-Va_BC) + IQ_BC*PetscCosScalar(Va_C-Va_BC)); /*dQ_C/dVang_C*/

        /* dP/dVmag */
        values[1]  = - (Vm_B/(PetscSqrtScalar(3.0)*Vm_AB))*(IP_AB*PetscCosScalar(Va_A+Va_B-2*Va_AB) - IQ_AB*PetscSinScalar(Va_A+Va_B-2*Va_AB))\
                     + (Vm_A/(PetscSqrtScalar(3.0)*Vm_AB*Vm_AB))*(Vm_A-Vm_B*PetscCosScalar(Va_A-Va_B))*(IP_AB*PetscCosScalar(Va_A-Va_AB) - IQ_AB*PetscSinScalar(Va_A-Va_AB))\
                     - (Vm_C/(PetscSqrtScalar(3.0)*Vm_CA))*(IP_CA*PetscCosScalar(Va_A+Va_C-2*Va_CA) - IQ_CA*PetscSinScalar(Va_A+Va_C-2*Va_CA))\
                     - (Vm_A/(PetscSqrtScalar(3.0)*Vm_CA*Vm_CA))*(Vm_A-Vm_C*PetscCosScalar(Va_C-Va_A))*(IP_CA*PetscCosScalar(Va_A-Va_CA) - IQ_CA*PetscSinScalar(Va_A-Va_CA)); /*dP_A/dVmag_A*/
        values[3]  = + (Vm_A/(PetscSqrtScalar(3.0)*Vm_AB))*(IP_AB*PetscCosScalar(Va_A+Va_B-2*Va_AB) - IQ_AB*PetscSinScalar(Va_A+Va_B-2*Va_AB))\
                     + (Vm_A/(PetscSqrtScalar(3.0)*Vm_AB*Vm_AB))*(Vm_B-Vm_A*PetscCosScalar(Va_A-Va_B))*(IP_AB*PetscCosScalar(Va_A-Va_AB) - IQ_AB*PetscSinScalar(Va_A-Va_AB)); /*dP_A/dVmag_B*/
        values[5]  = + (Vm_A/(PetscSqrtScalar(3.0)*Vm_CA))*(IP_CA*PetscCosScalar(Va_C+Va_A-2*Va_CA) - IQ_CA*PetscSinScalar(Va_C+Va_A-2*Va_CA))\
                     - (Vm_A/(PetscSqrtScalar(3.0)*Vm_CA*Vm_CA))*(Vm_C-Vm_A*PetscCosScalar(Va_C-Va_A))*(IP_CA*PetscCosScalar(Va_A-Va_CA) - IQ_CA*PetscSinScalar(Va_A-Va_CA)); /*dP_A/dVmag_C*/

        values[13] = + (Vm_B/(PetscSqrtScalar(3.0)*Vm_AB))*(IP_AB*PetscCosScalar(Va_A+Va_B-2*Va_AB) - IQ_AB*PetscSinScalar(Va_A+Va_B-2*Va_AB))\
                     - (Vm_B/(PetscSqrtScalar(3.0)*Vm_AB*Vm_AB))*(Vm_A-Vm_B*PetscCosScalar(Va_A-Va_B))*(IP_AB*PetscCosScalar(Va_B-Va_AB) - IQ_AB*PetscSinScalar(Va_B-Va_AB)); /*dP_B/dVmag_A*/
        values[15] = - (Vm_C/(PetscSqrtScalar(3.0)*Vm_BC))*(IP_BC*PetscCosScalar(Va_B+Va_C-2*Va_BC) - IQ_BC*PetscSinScalar(Va_B+Va_C-2*Va_BC))\
                     + (Vm_B/(PetscSqrtScalar(3.0)*Vm_BC*Vm_BC))*(Vm_B-Vm_C*PetscCosScalar(Va_B-Va_C))*(IP_BC*PetscCosScalar(Va_B-Va_BC) - IQ_BC*PetscSinScalar(Va_B-Va_BC))\
                     - (Vm_A/(PetscSqrtScalar(3.0)*Vm_AB))*(IP_AB*PetscCosScalar(Va_A+Va_B-2*Va_AB) - IQ_AB*PetscSinScalar(Va_A+Va_B-2*Va_AB))\
                     - (Vm_B/(PetscSqrtScalar(3.0)*Vm_AB*Vm_AB))*(Vm_B-Vm_A*PetscCosScalar(Va_A-Va_B))*(IP_AB*PetscCosScalar(Va_B-Va_AB) - IQ_AB*PetscSinScalar(Va_B-Va_AB)); /*dP_B/dVmag_B*/
        values[17] = + (Vm_B/(PetscSqrtScalar(3.0)*Vm_BC))*(IP_BC*PetscCosScalar(Va_B+Va_C-2*Va_BC) - IQ_BC*PetscSinScalar(Va_B+Va_C-2*Va_BC))\
                     + (Vm_B/(PetscSqrtScalar(3.0)*Vm_BC*Vm_BC))*(Vm_C-Vm_B*PetscCosScalar(Va_B-Va_C))*(IP_BC*PetscCosScalar(Va_B-Va_BC) - IQ_BC*PetscSinScalar(Va_B-Va_BC)); /*dP_B/dVmag_C*/

        values[25] = + (Vm_C/(PetscSqrtScalar(3.0)*Vm_CA))*(IP_CA*PetscCosScalar(Va_C+Va_A-2*Va_CA) - IQ_CA*PetscSinScalar(Va_C+Va_A-2*Va_CA))\
                     + (Vm_C/(PetscSqrtScalar(3.0)*Vm_CA*Vm_CA))*(Vm_A-Vm_C*PetscCosScalar(Va_C-Va_A))*(IP_CA*PetscCosScalar(Va_C-Va_CA) - IQ_CA*PetscSinScalar(Va_C-Va_CA)); /*dP_C/dVmag_A*/
        values[27] = + (Vm_C/(PetscSqrtScalar(3.0)*Vm_BC))*(IP_BC*PetscCosScalar(Va_B+Va_C-2*Va_BC) - IQ_BC*PetscSinScalar(Va_B+Va_C-2*Va_BC))\
                     - (Vm_C/(PetscSqrtScalar(3.0)*Vm_BC*Vm_BC))*(Vm_B-Vm_C*PetscCosScalar(Va_B-Va_C))*(IP_BC*PetscCosScalar(Va_C-Va_BC) - IQ_BC*PetscSinScalar(Va_C-Va_BC)); /*dP_C/dVmag_B*/
        values[29] = - (Vm_A/(PetscSqrtScalar(3.0)*Vm_CA))*(IP_CA*PetscCosScalar(Va_C+Va_A-2*Va_CA) - IQ_CA*PetscSinScalar(Va_C+Va_A-2*Va_CA))\
                     + (Vm_C/(PetscSqrtScalar(3.0)*Vm_CA*Vm_CA))*(Vm_C-Vm_A*PetscCosScalar(Va_C-Va_A))*(IP_CA*PetscCosScalar(Va_C-Va_CA) - IQ_CA*PetscSinScalar(Va_C-Va_CA))\
                     - (Vm_B/(PetscSqrtScalar(3.0)*Vm_BC))*(IP_BC*PetscCosScalar(Va_B+Va_C-2*Va_BC) - IQ_BC*PetscSinScalar(Va_B+Va_C-2*Va_BC))\
                     - (Vm_C/(PetscSqrtScalar(3.0)*Vm_BC*Vm_BC))*(Vm_C-Vm_B*PetscCosScalar(Va_B-Va_C))*(IP_BC*PetscCosScalar(Va_C-Va_BC) - IQ_BC*PetscSinScalar(Va_C-Va_BC)); /*dP_C/dVmag_C*/

        /* dQ/dVmag */
        values[7]  = - (Vm_B/(PetscSqrtScalar(3.0)*Vm_AB))*(IP_AB*PetscSinScalar(Va_A+Va_B-2*Va_AB) + IQ_AB*PetscCosScalar(Va_A+Va_B-2*Va_AB))\
                     + (Vm_A/(PetscSqrtScalar(3.0)*Vm_AB*Vm_AB))*(Vm_A-Vm_B*PetscCosScalar(Va_A-Va_B))*(IP_AB*PetscSinScalar(Va_A-Va_AB) + IQ_AB*PetscCosScalar(Va_A-Va_AB))\
                     - (Vm_C/(PetscSqrtScalar(3.0)*Vm_CA))*(IP_CA*PetscSinScalar(Va_A+Va_C-2*Va_CA) + IQ_CA*PetscCosScalar(Va_A+Va_C-2*Va_CA))\
                     - (Vm_A/(PetscSqrtScalar(3.0)*Vm_CA*Vm_CA))*(Vm_A-Vm_C*PetscCosScalar(Va_C-Va_A))*(IP_CA*PetscSinScalar(Va_A-Va_CA) + IQ_CA*PetscCosScalar(Va_A-Va_CA)); /*dQ_A/dVmag_A*/
        values[9]  = + (Vm_A/(PetscSqrtScalar(3.0)*Vm_AB))*(IP_AB*PetscSinScalar(Va_A+Va_B-2*Va_AB) + IQ_AB*PetscCosScalar(Va_A+Va_B-2*Va_AB))\
                     + (Vm_A/(PetscSqrtScalar(3.0)*Vm_AB*Vm_AB))*(Vm_B-Vm_A*PetscCosScalar(Va_A-Va_B))*(IP_AB*PetscSinScalar(Va_A-Va_AB) + IQ_AB*PetscCosScalar(Va_A-Va_AB)); /*dQ_A/dVmag_B*/
        values[11] = + (Vm_A/(PetscSqrtScalar(3.0)*Vm_CA))*(IP_CA*PetscSinScalar(Va_C+Va_A-2*Va_CA) + IQ_CA*PetscCosScalar(Va_C+Va_A-2*Va_CA))\
                     - (Vm_A/(PetscSqrtScalar(3.0)*Vm_CA*Vm_CA))*(Vm_C-Vm_A*PetscCosScalar(Va_C-Va_A))*(IP_CA*PetscSinScalar(Va_A - Va_CA) + IQ_CA*PetscCosScalar(Va_A - Va_CA)); /*dQ_A/dVmag_C*/

        values[19] = + (Vm_B/(PetscSqrtScalar(3.0)*Vm_AB))*(IP_AB*PetscSinScalar(Va_A+Va_B-2*Va_AB) + IQ_AB*PetscCosScalar(Va_A+Va_B-2*Va_AB))\
                     - (Vm_B/(PetscSqrtScalar(3.0)*Vm_AB*Vm_AB))*(Vm_A-Vm_B*PetscCosScalar(Va_A-Va_B))*(IP_AB*PetscSinScalar(Va_B-Va_AB) + IQ_AB*PetscCosScalar(Va_B-Va_AB)); /*dQ_B/dVmag_A*/
        values[21] = - (Vm_C/(PetscSqrtScalar(3.0)*Vm_BC))*(IP_BC*PetscSinScalar(Va_B+Va_C-2*Va_BC) + IQ_BC*PetscCosScalar(Va_B+Va_C-2*Va_BC))\
                     + (Vm_B/(PetscSqrtScalar(3.0)*Vm_BC*Vm_BC))*(Vm_B-Vm_C*PetscCosScalar(Va_B-Va_C))*(IP_BC*PetscSinScalar(Va_B-Va_BC) + IQ_BC*PetscCosScalar(Va_B-Va_BC))\
                     - (Vm_A/(PetscSqrtScalar(3.0)*Vm_AB))*(IP_AB*PetscSinScalar(Va_A+Va_B-2*Va_AB) + IQ_AB*PetscCosScalar(Va_A+Va_B-2*Va_AB))\
                     - (Vm_B/(PetscSqrtScalar(3.0)*Vm_AB*Vm_AB))*(Vm_B-Vm_A*PetscCosScalar(Va_A-Va_B))*(IP_AB*PetscSinScalar(Va_B-Va_AB) + IQ_AB*PetscCosScalar(Va_B-Va_AB)); /*dQ_B/dVmag_B*/
        values[23] = + (Vm_B/(PetscSqrtScalar(3.0)*Vm_BC))*(IP_BC*PetscSinScalar(Va_B+Va_C-2*Va_BC) + IQ_BC*PetscCosScalar(Va_B+Va_C-2*Va_BC))\
                     + (Vm_B/(PetscSqrtScalar(3.0)*Vm_BC*Vm_BC))*(Vm_C-Vm_B*PetscCosScalar(Va_B-Va_C))*(IP_BC*PetscSinScalar(Va_B-Va_BC) + IQ_BC*PetscCosScalar(Va_B-Va_BC)); /*dQ_B/dVmag_C*/

        values[31] = + (Vm_C/(PetscSqrtScalar(3.0)*Vm_CA))*(IP_CA*PetscSinScalar(Va_C+Va_A-2*Va_CA) + IQ_CA*PetscCosScalar(Va_C+Va_A-2*Va_CA))\
                     + (Vm_C/(PetscSqrtScalar(3.0)*Vm_CA*Vm_CA))*(Vm_A-Vm_C*PetscCosScalar(Va_C-Va_A))*(IP_CA*PetscSinScalar(Va_C-Va_CA) + IQ_CA*PetscCosScalar(Va_C-Va_CA)); /*dQ_C/dVmag_A*/
        values[33] = + (Vm_C/(PetscSqrtScalar(3.0)*Vm_BC))*(IP_BC*PetscSinScalar(Va_B+Va_C-2*Va_BC) + IQ_BC*PetscCosScalar(Va_B+Va_C-2*Va_BC))\
                     - (Vm_C/(PetscSqrtScalar(3.0)*Vm_BC*Vm_BC))*(Vm_B-Vm_C*PetscCosScalar(Va_B-Va_C))*(IP_BC*PetscSinScalar(Va_C-Va_BC) + IQ_BC*PetscCosScalar(Va_C-Va_BC)); /*dQ_C/dVmag_B*/
        values[35] = - (Vm_A/(PetscSqrtScalar(3.0)*Vm_CA))*(IP_CA*PetscSinScalar(Va_C+Va_A-2*Va_CA) + IQ_CA*PetscCosScalar(Va_C+Va_A-2*Va_CA))\
                     + (Vm_C/(PetscSqrtScalar(3.0)*Vm_CA*Vm_CA))*(Vm_C-Vm_A*PetscCosScalar(Va_C-Va_A))*(IP_CA*PetscSinScalar(Va_C-Va_CA) + IQ_CA*PetscCosScalar(Va_C-Va_CA))\
                     - (Vm_B/(PetscSqrtScalar(3.0)*Vm_BC))*(IP_BC*PetscSinScalar(Va_B+Va_C-2*Va_BC) + IQ_BC*PetscCosScalar(Va_B+Va_C-2*Va_BC))\
                     - (Vm_C/(PetscSqrtScalar(3.0)*Vm_BC*Vm_BC))*(Vm_C-Vm_B*PetscCosScalar(Va_B-Va_C))*(IP_BC*PetscSinScalar(Va_C-Va_BC) + IQ_BC*PetscCosScalar(Va_C-Va_BC)); /*dQ_C/dVmag_C*/

        ierr = MatSetValues(J,6,row,6,col,values,ADD_VALUES);CHKERRQ(ierr);

      } /* End of conditional for constant-current delta loads */

    } /* End of Jacobian modification because of delta loads */

  } /* End of vertices for loop */
  ierr = VecRestoreArrayRead(localX,&xarr);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(networkdm,&localX);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
/*
  PetscInt p,q,rows,cols;
  PetscViewer viewer;
  PetscScalar val;
  ierr = MatGetSize(J,&rows,&cols);CHKERRQ(ierr);
  ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"JacobianFile",&viewer);CHKERRQ(ierr);
  for(p=0; p<rows; p++){
    for(q=0; q<cols; q++){
      ierr = MatGetValues(J,1,&p,1,&q,&val);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(viewer,"%.4lf\t",val);CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
*/
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SetInitialValues"
PetscErrorCode SetInitialValues(DM networkdm,Vec X,void* appctx)
{
  PetscErrorCode ierr;
  VERTEXDATA     bus;
  PetscInt       v, vStart, vEnd, offset;
  Vec            localX;
  PetscScalar    *xarr;
  PetscInt       key,numComps,j,offsetd;
  DMNetworkComponentGenericDataType *arr;

  PetscFunctionBegin;
  ierr = DMNetworkGetVertexRange(networkdm,&vStart, &vEnd);CHKERRQ(ierr);

  ierr = DMGetLocalVector(networkdm,&localX);CHKERRQ(ierr);

  ierr = VecSet(X,0.0);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  ierr = VecGetArray(localX,&xarr);CHKERRQ(ierr);
  ierr = DMNetworkGetComponentDataArray(networkdm,&arr);CHKERRQ(ierr);
  for (v = vStart; v < vEnd; v++) {
    ierr = DMNetworkGetVariableOffset(networkdm,v,&offset);CHKERRQ(ierr);
    ierr = DMNetworkGetNumComponents(networkdm,v,&numComps);CHKERRQ(ierr);
    for (j=0; j < numComps; j++) {
      ierr = DMNetworkGetComponentTypeOffset(networkdm,v,j,&key,&offsetd);CHKERRQ(ierr);
      if (key == 1) {
	       bus = (VERTEXDATA)(arr+offsetd);
	       xarr[offset] = bus->va*PETSC_PI/180.0;
	       xarr[offset+1] = bus->vm;
      }
    }
  }
  ierr = VecRestoreArray(localX,&xarr);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(networkdm,localX,ADD_VALUES,X);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(networkdm,localX,ADD_VALUES,X);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(networkdm,&localX);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PrintDPFOutput"
PetscErrorCode PrintDPFOutput(Vec X)
{
  PetscViewer    viewer;
  PetscErrorCode ierr;
  PetscInt       X_size, i;
  const PetscScalar *xarr;

  PetscFunctionBegin;

  ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"DPF_output",&viewer);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"Bus data:\n");CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"VA_mag\t\t\tVA_ang\t\t\tVB_mag\t\t\tVB_ang\t\t\tVC_mag\t\t\tVC_ang\n");CHKERRQ(ierr);

  ierr = VecGetArrayRead(X,&xarr);CHKERRQ(ierr);
  ierr = VecGetSize(X,&X_size);CHKERRQ(ierr);

  for (i=0; i<X_size/6; i++){
    ierr = PetscViewerASCIIPrintf(viewer,"%lf\t\t%lf\t\t%lf\t\t%lf\t\t%lf\t\t%lf\n",\
      xarr[6*i+1],xarr[6*i]*180.0/PETSC_PI,xarr[6*i+3],xarr[6*i+2]*180.0/PETSC_PI,\
      xarr[6*i+5],xarr[6*i+4]*180.0/PETSC_PI);CHKERRQ(ierr);
  }

  ierr = VecRestoreArrayRead(X,&xarr);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char ** argv)
{
  PetscErrorCode ierr;
  char           pfdata_file[PETSC_MAX_PATH_LEN]="datafiles/distCase_4Dyn1.m";
  PFDATA         *pfdata;
  PetscInt       numEdges=0,numVertices=0;
  int            *edges = NULL;
  PetscInt       i;
  DM             networkdm;
  PetscInt       componentkey[5];
  UserCtx        User;
  PetscLogStage  stage1,stage2;
  PetscInt       eStart, eEnd, vStart, vEnd,j;
  PetscInt       genj,loadj,dloadj;
  Vec            X,F,X_local;
  Mat            J;
  SNES           snes;
  VecScatter     scatter_ctx;
  PetscScalar    t1, t2;

  PetscInitialize(&argc,&argv,"pfoptions",help);
  ierr = PetscTime(&t1);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

  /* Create an empty network object */
  ierr = DMNetworkCreate(PETSC_COMM_WORLD,&networkdm);CHKERRQ(ierr);
  /* Register the components in the network */
  ierr = DMNetworkRegisterComponent(networkdm,"branchstruct",sizeof(struct _p_EDGEDATA),&componentkey[0]);CHKERRQ(ierr);
  ierr = DMNetworkRegisterComponent(networkdm,"busstruct",sizeof(struct _p_VERTEXDATA),&componentkey[1]);CHKERRQ(ierr);
  ierr = DMNetworkRegisterComponent(networkdm,"genstruct",sizeof(struct _p_GEN),&componentkey[2]);CHKERRQ(ierr);
  ierr = DMNetworkRegisterComponent(networkdm,"wyeloadstruct",sizeof(struct _p_WYELOAD),&componentkey[3]);CHKERRQ(ierr);
  ierr = DMNetworkRegisterComponent(networkdm,"deltaloadstruct",sizeof(struct _p_DELTALOAD),&componentkey[4]);CHKERRQ(ierr);


/******************************************************** Stage 1: data lecture ********************************************************/
  ierr = PetscLogStageRegister("Read Data",&stage1);CHKERRQ(ierr);
  PetscLogStagePush(stage1);
  /* READ THE DATA */
  if (!rank) {
    /*    READ DATA */
    /* Only rank 0 reads the data */
    ierr = PetscOptionsGetString(NULL,NULL,"-pfdata",pfdata_file,PETSC_MAX_PATH_LEN-1,NULL);CHKERRQ(ierr);
    ierr = PetscNew(&pfdata);CHKERRQ(ierr);
    ierr = DPFReadMatPowerData(pfdata,pfdata_file);CHKERRQ(ierr);

    ierr = PrintDPF(pfdata);CHKERRQ(ierr);

    User.Sbase = pfdata->sbase;
    numEdges = pfdata->nbranch;
    numVertices = pfdata->nbus;

    ierr = PetscMalloc1(2*numEdges,&edges);CHKERRQ(ierr);
    ierr = GetListofEdges(pfdata->nbranch,pfdata->branch,edges);CHKERRQ(ierr);
  }
  PetscLogStagePop();
/************************************************************* End of stage 1 *************************************************************/

  ierr = MPI_Barrier(PETSC_COMM_WORLD);CHKERRQ(ierr);

/******************************************************** Stage 2: network creation ********************************************************/
  ierr = PetscLogStageRegister("Create network",&stage2);CHKERRQ(ierr);
  PetscLogStagePush(stage2);
  /* Set number of nodes/edges */
  ierr = DMNetworkSetSizes(networkdm,numVertices,numEdges,PETSC_DETERMINE,PETSC_DETERMINE);CHKERRQ(ierr);
  /* Add edge connectivity */
  ierr = DMNetworkSetEdgeList(networkdm,edges);CHKERRQ(ierr);
  /* Set up the network layout */
  ierr = DMNetworkLayoutSetUp(networkdm);CHKERRQ(ierr);

  if (!rank) {
    ierr = PetscFree(edges);CHKERRQ(ierr);
  }
  /* Add network components */

  genj=0; loadj=0; dloadj=0;

  ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);
  for (i = eStart; i < eEnd; i++) {
    ierr = DMNetworkAddComponent(networkdm,i,componentkey[0],&pfdata->branch[i-eStart]);CHKERRQ(ierr);
  }

  ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr);
  for (i = vStart; i < vEnd; i++){
    ierr = DMNetworkAddComponent(networkdm,i,componentkey[1],&pfdata->bus[i-vStart]);CHKERRQ(ierr);
    if (pfdata->bus[i-vStart].ngen) {
      for (j = 0; j < pfdata->bus[i-vStart].ngen; j++) {
	       ierr = DMNetworkAddComponent(networkdm,i,componentkey[2],&pfdata->gen[genj++]);CHKERRQ(ierr);
      }
    }
    if (pfdata->bus[i-vStart].nload) {
      for (j=0; j < pfdata->bus[i-vStart].nload; j++) {
	       ierr = DMNetworkAddComponent(networkdm,i,componentkey[3],&pfdata->wye[loadj++]);CHKERRQ(ierr);
      }
    }
    if (pfdata->bus[i-vStart].ndload) {
      for (j=0; j < pfdata->bus[i-vStart].ndload; j++) {
	       ierr = DMNetworkAddComponent(networkdm,i,componentkey[4],&pfdata->delta[dloadj++]);CHKERRQ(ierr);
      }
    }
    /* Add number of variables */
    ierr = DMNetworkAddNumVariables(networkdm,i,2);CHKERRQ(ierr);
  }

  /* Set up DM for use */
  ierr = DMSetUp(networkdm);CHKERRQ(ierr);

  if (!rank) {
    ierr = PetscFree(pfdata->bus);CHKERRQ(ierr);
    ierr = PetscFree(pfdata->gen);CHKERRQ(ierr);
    ierr = PetscFree(pfdata->branch);CHKERRQ(ierr);
    ierr = PetscFree(pfdata->wye);CHKERRQ(ierr);
    ierr = PetscFree(pfdata->delta);CHKERRQ(ierr);
    ierr = PetscFree(pfdata);CHKERRQ(ierr);
  }

  PetscLogStagePop();
/************************************************************* End of stage 2 *************************************************************/
/*
  ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);
  ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_SELF,"%d: %d,%d,%d,%d\n",rank,eStart,eEnd,vStart,vEnd);CHKERRQ(ierr);
*/
  /* Broadcast Sbase to all processors */
  ierr = MPI_Bcast(&User.Sbase,1,MPIU_SCALAR,0,PETSC_COMM_WORLD);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(networkdm,&X);CHKERRQ(ierr);
  ierr = VecDuplicate(X,&F);CHKERRQ(ierr);

  ierr = DMCreateMatrix(networkdm,&J);CHKERRQ(ierr);
  ierr = MatSetOption(J,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE);CHKERRQ(ierr);

  ierr = SetInitialValues(networkdm,X,&User);CHKERRQ(ierr);

  /* HOOK UP SOLVER */
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
  ierr = SNESSetDM(snes,networkdm);CHKERRQ(ierr);
  ierr = SNESSetFunction(snes,F,FormFunction,&User);CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,J,J,FormJacobian,&User);CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  ierr = SNESSolve(snes,NULL,X);CHKERRQ(ierr);

  ierr = VecScatterCreateToZero(X,&scatter_ctx,&X_local);CHKERRQ(ierr);
  ierr = VecScatterBegin(scatter_ctx,X,X_local,INSERT_VALUES,SCATTER_FORWARD);
  ierr = VecScatterEnd(scatter_ctx,X,X_local,INSERT_VALUES,SCATTER_FORWARD);
  ierr = VecScatterDestroy(&scatter_ctx);

  if (!rank) {
    ierr = PrintDPFOutput(X_local);CHKERRQ(ierr);
  }
  ierr = VecDestroy(&X);CHKERRQ(ierr);
  ierr = VecDestroy(&X_local);CHKERRQ(ierr);
  ierr = VecDestroy(&F);CHKERRQ(ierr);
  ierr = MatDestroy(&J);CHKERRQ(ierr);
  ierr = SNESDestroy(&snes);CHKERRQ(ierr);
  ierr = DMDestroy(&networkdm);CHKERRQ(ierr);

  ierr = PetscTime(&t2);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Total time elapsed: %lf seconds\n",t2-t1);CHKERRQ(ierr);

  PetscFinalize();
  return 0;
}
